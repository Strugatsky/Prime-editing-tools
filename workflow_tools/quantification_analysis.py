import argparse
import sqlite3
import pandas as pd
import re
import sys
import hashlib
import uuid

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process prime editing data and update SQLite database.')
    parser.add_argument('tsv_file', type=str, help='Path to the TSV file containing editing data')
    parser.add_argument('db_file', type=str, help='Path to the SQLite database')
    return parser.parse_args()

def extract_batch_info(batch_name):
    # PE, P  R
    pattern1 = r'[a-zA-Z]+(PE\w+)_P(\d+)_R(\d+)_[R/r]ep(\d+)'
    # Match pattern with PE type, R first then P
    pattern2 = r'[a-zA-Z]+(PE\w+)_R(\d+)_P(\d+)_[R/r]ep(\d+)'
    # Match pattern without PE type, P first then R
    pattern3 = r'[a-zA-Z]+_P(\d+)_R(\d+)_[R/r]ep(\d+)'
    # Match pattern without PE type, R first then P
    pattern4 = r'[a-zA-Z]+_R(\d+)_P(\d+)_[R/r]ep(\d+)'
    # Match pattern with drug suffix (no rep)
    pattern_drug = r'[a-zA-Z_]+(PE\w+)_P(\d+)R(\d+)_([a-zA-Z]+)$'
    # _drug_rep
    pattern_drug_rep = r'[a-zA-Z_]+(PE\w+)_P(\d+)R(\d+)_([a-zA-Z]+)_[R/r]ep(\d+)'

    
        
    # Try other patterns (no drug)
    match1 = re.match(pattern1, batch_name)
    match2 = re.match(pattern2, batch_name)
    match3 = re.match(pattern3, batch_name)
    match4 = re.match(pattern4, batch_name)
    match_drug_rep = re.match(pattern_drug_rep, batch_name)
    match_drug = re.match(pattern_drug, batch_name)
    
    if match1:
        prime_editor = match1.group(1)
        pbs = int(match1.group(2))
        rtt = int(match1.group(3))
        replicate = int(match1.group(4))
        return prime_editor, pbs, rtt, replicate, None
    elif match2:
        prime_editor = match2.group(1)
        rtt = int(match2.group(2))
        pbs = int(match2.group(3))
        replicate = int(match2.group(4))
        return prime_editor, pbs, rtt, replicate, None
    elif match3:
        pbs = int(match3.group(1))
        rtt = int(match3.group(2))
        replicate = int(match3.group(3))
        return None, pbs, rtt, replicate, None
    elif match4:
        rtt = int(match4.group(1))
        pbs = int(match4.group(2))
        replicate = int(match4.group(3))
        return None, pbs, rtt, replicate, None
    elif match_drug_rep:
        prime_editor = match_drug_rep.group(1)
        pbs = int(match_drug_rep.group(2))
        rtt = int(match_drug_rep.group(3))
        replicate = int(match_drug_rep.group(4))
        drug_code = match_drug_rep.group(5)
        return prime_editor, pbs, rtt, replicate, drug_code
    elif match_drug:
        prime_editor = match_drug.group(1)
        pbs = int(match_drug.group(2))
        rtt = int(match_drug.group(3))
        drug_code = match_drug.group(4)
        replicate = 1
        return prime_editor, pbs, rtt, replicate, drug_code
    else:
        print(f"Warning: Could not parse batch name: {batch_name}")
        return None, None, None, None, None

def process_data(tsv_file):
    """Process the TSV file to calculate editing efficiencies for each sample."""
    # Read TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Group by batch name (3 rows per batch)
    batches = {}
    unique_batches = df['Batch'].unique()
    
    # Dictionary to track batches with missing prime editor
    missing_pe_batches = {}
    
    for batch_name in unique_batches:
        batch_data = df[df['Batch'] == batch_name]
        
        if len(batch_data) != 3:
            print(f"Warning: Batch {batch_name} does not have exactly 3 rows. Skipping.")
            continue
        
        # Calculate total reads for this batch
        total_unmodified = batch_data['Unmodified'].sum()
        total_modified = batch_data['Modified'].sum()
        total_discarded = batch_data['Discarded'].sum()
        total_reads = total_unmodified + total_modified + total_discarded
        
        # Get prime-edited row
        prime_edited_row = batch_data[batch_data['Amplicon'] == 'Prime-edited']
        if len(prime_edited_row) != 1:
            print(f"Warning: Batch {batch_name} does not have exactly one 'Prime-edited' row. Skipping.")
            continue
            
        # Get scaffold-incorporated row
        scaffold_row = batch_data[batch_data['Amplicon'] == 'Scaffold-incorporated']
        if len(scaffold_row) != 1:
            print(f"Warning: Batch {batch_name} does not have exactly one 'Scaffold-incorporated' row. Skipping.")
            continue
        
        # Calculate percentages
        correct_edits = (prime_edited_row['Unmodified'].values[0] / total_reads) * 100
        incorrect_edits = (prime_edited_row['Modified'].values[0] / total_reads) * 100
        scaffold_incorporated = (scaffold_row['Modified'].values[0] / total_reads) * 100
        
        # Extract batch info
        prime_editor, pbs, rtt, replicate, drug_code = extract_batch_info(batch_name)
        
        if None in (pbs, rtt, replicate):
            print(f"Warning: Could not extract PBS, RTT, or replicate from batch {batch_name}. Skipping.")
            continue
        
        # Check for control samples
        if drug_code and drug_code.lower() == 'ctrl':
            drug_code = None
        
        # If prime editor is missing, track it for later user input
        if prime_editor is None:
            missing_pe_batches[batch_name] = {
                'pbs': pbs,
                'rtt': rtt,
                'replicate': replicate,
                'correct_edits': correct_edits,
                'incorrect_edits': incorrect_edits,
                'scaffold_incorporated': scaffold_incorporated,
                'drug_code': drug_code
            }
        else:
            batches[batch_name] = {
                'prime_editor': prime_editor,
                'pbs': pbs,
                'rtt': rtt,
                'replicate': replicate,
                'correct_edits': correct_edits,
                'incorrect_edits': incorrect_edits,
                'scaffold_incorporated': scaffold_incorporated,
                'drug_code': drug_code
            }
    
    # Handle batches with missing prime editor
    if missing_pe_batches:
        print(f"\nFound {len(missing_pe_batches)} batches with missing prime editor information.")
        response = input("Enter prime editor name for these batches (e.g., PE2, PEMax), or enter different for individual prompts: ")
        
        if response.lower() == 'different':
            # Ask for each batch individually
            for batch_name, data in missing_pe_batches.items():
                prime_editor = input(f"Enter prime editor name for batch '{batch_name}': ")
                data['prime_editor'] = prime_editor
                batches[batch_name] = data
        else:
            # Use the same prime editor for all missing batches
            for batch_name, data in missing_pe_batches.items():
                data['prime_editor'] = response
                batches[batch_name] = data
    
    return batches

def setup_database(conn):
    """Set up the new tables in the database if they don't exist."""
    cursor = conn.cursor()
    
    # Create runs table if it doesn't exist
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS runs (
        id TEXT PRIMARY KEY,
        run_name TEXT NOT NULL,
        experiment_id TEXT NOT NULL,
        FOREIGN KEY (experiment_id) REFERENCES experiments(id)
    )
    ''')
    
    # Create drugs table if it doesn't exist
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS drugs (
        id TEXT PRIMARY KEY,
        name TEXT NOT NULL,
        description TEXT
    )
    ''')
    
    # Insert default "None" drug if it doesn't exist
    cursor.execute('''
    INSERT OR IGNORE INTO drugs (id, name, description)
    VALUES (?, ?, ?)
    ''', ('00000000-0000-0000-0000-000000000000', 'None', 'No drug used'))
    
    # Check if the drug_id column exists in data_points
    cursor.execute("PRAGMA table_info(data_points)")
    columns = cursor.fetchall()
    column_names = [column[1] for column in columns]
    
    if 'data_points' in get_table_names(conn) and 'drug_id' not in column_names:
        # Add drug_id column to existing data_points table
        cursor.execute('''
        ALTER TABLE data_points
        ADD COLUMN drug_id TEXT DEFAULT '00000000-0000-0000-0000-000000000000'
        REFERENCES drugs(id)
        ''')
        print("Added drug_id column to existing data_points table")
    
    # Create data_points table if it doesn't exist (now with drug_id)
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS data_points (
        id TEXT PRIMARY KEY,
        experiment_entry_id TEXT NOT NULL,
        correct_edits REAL NOT NULL,
        incorrect_edits REAL NOT NULL,
        scaffold_incorporated REAL NOT NULL,
        prime_editor TEXT NOT NULL,
        replicate INTEGER NOT NULL,
        run_id TEXT NOT NULL,
        drug_id TEXT DEFAULT '00000000-0000-0000-0000-000000000000',
        FOREIGN KEY (experiment_entry_id) REFERENCES experiment_entries(id),
        FOREIGN KEY (run_id) REFERENCES runs(id),
        FOREIGN KEY (drug_id) REFERENCES drugs(id)
    )
    ''')
    
    conn.commit()

def get_table_names(conn):
    """Get a list of all table names in the database."""
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    return [table[0] for table in cursor.fetchall()]

def get_experiment_entry_id(conn, experiment_id, pbs, rtt):
    """Find the experiment entry ID based on experiment_id, PBS, and RTT values."""
    cursor = conn.cursor()
    cursor.execute('''
    SELECT id FROM experiment_entries 
    WHERE experiment_id = ? AND pbs = ? AND rtt = ?
    ''', (experiment_id, pbs, rtt))
    
    result = cursor.fetchone()
    if result:
        return result[0]
    return None

def get_available_experiments(conn):
    """Get a list of available experiments from the database."""
    cursor = conn.cursor()
    cursor.execute('''
    SELECT id, name, date, variant FROM experiments
    ''')
    
    return cursor.fetchall()

def get_or_create_drug(conn, drug_code):
    """Get an existing drug or create a new one based on drug code."""
    if drug_code is None:
        # Return the default "None" drug
        return '00000000-0000-0000-0000-000000000000'
    
    cursor = conn.cursor()
    
    # Check if we've seen this drug code before in this run
    if drug_code in get_or_create_drug.drug_cache:
        return get_or_create_drug.drug_cache[drug_code]
    
    # Get all available drugs
    cursor.execute('''
    SELECT id, name FROM drugs
    ''')
    drugs = cursor.fetchall()
    
    print(f"\nFound drug code '{drug_code}' in batch name.")
    print("Available drugs:")
    print("0. Create new drug")
    for i, (drug_id, drug_name) in enumerate(drugs, 1):
        print(f"{i}. {drug_name}")
    
    while True:
        try:
            selection = int(input("\nSelect a drug or create new (0-{}): ".format(len(drugs))))
            if 0 <= selection <= len(drugs):
                break
            print("Invalid selection. Please try again.")
        except ValueError:
            print("Please enter a number.")
    
    if selection == 0:
        # Create new drug
        drug_name = input(f"Enter name for drug '{drug_code}': ")
        drug_description = input("Enter description (optional): ")
        drug_id = str(uuid.uuid4())
        
        cursor.execute('''
        INSERT INTO drugs (id, name, description)
        VALUES (?, ?, ?)
        ''', (drug_id, drug_name, drug_description))
        conn.commit()
        
    else:
        # Use existing drug
        drug_id = drugs[selection-1][0]
    
    # Cache this drug code for future use in this run
    get_or_create_drug.drug_cache[drug_code] = drug_id
    return drug_id

# Initialize drug cache
get_or_create_drug.drug_cache = {}

def insert_run(conn, run_name, experiment_id):
    """Insert a new run into the runs table."""
    cursor = conn.cursor()
    run_id = str(uuid.uuid4())
    
    cursor.execute('''
    INSERT INTO runs (id, run_name, experiment_id)
    VALUES (?, ?, ?)
    ''', (run_id, run_name, experiment_id))
    
    conn.commit()
    return run_id

def insert_data_point(conn, experiment_entry_id, correct_edits, incorrect_edits, 
                     scaffold_incorporated, prime_editor, replicate, run_id, drug_id):
    """Insert a data point into the data_points table."""
    cursor = conn.cursor()
    data_point_id = str(uuid.uuid4())
    
    cursor.execute('''
    INSERT INTO data_points (id, experiment_entry_id, correct_edits, incorrect_edits,
                           scaffold_incorporated, prime_editor, replicate, run_id, drug_id)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (data_point_id, experiment_entry_id, correct_edits, incorrect_edits,
         scaffold_incorporated, prime_editor, replicate, run_id, drug_id))
    
    conn.commit()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Process the TSV file
    try:
        batches = process_data(args.tsv_file)
    except Exception as e:
        print(f"Error processing TSV file: {e}")
        return 1
    
    # Connect to the database
    try:
        conn = sqlite3.connect(args.db_file)
        setup_database(conn)
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return 1
    
    # Display available experiments
    experiments = get_available_experiments(conn)
    if not experiments:
        print("No experiments found in the database.")
        conn.close()
        return 1
    
    print("Available experiments:")
    for i, (exp_id, name, date, variant) in enumerate(experiments, 1):
        print(f"{i}. {name} ({date}) - {variant}")
    
    # Get user input for experiment selection
    while True:
        try:
            selection = int(input("\nSelect an experiment (1-{}): ".format(len(experiments))))
            if 1 <= selection <= len(experiments):
                break
            print("Invalid selection. Please try again.")
        except ValueError:
            print("Please enter a number.")
    
    selected_experiment_id = experiments[selection-1][0]
    selected_experiment_name = experiments[selection-1][1]
    
    # Get run name from user
    run_name = input(f"\nEnter a run name for experiment '{selected_experiment_name}': ")
    
    # Insert run
    run_id = insert_run(conn, run_name, selected_experiment_id)
    
    # Reset drug cache for this run
    get_or_create_drug.drug_cache = {}
    
    # Process batch data and insert into database
    data_points_added = 0
    for batch_name, batch_data in batches.items():
        # Get experiment entry ID
        experiment_entry_id = get_experiment_entry_id(
            conn, 
            selected_experiment_id, 
            batch_data['pbs'], 
            batch_data['rtt']
        )
        
        if not experiment_entry_id:
            print(f"Warning: No experiment entry found for PBS={batch_data['pbs']}, RTT={batch_data['rtt']}. Skipping batch {batch_name}.")
            continue
        
        # Handle drug if present
        drug_id = get_or_create_drug(conn, batch_data.get('drug_code'))
        
        # Insert data point
        insert_data_point(
            conn,
            experiment_entry_id,
            batch_data['correct_edits'],
            batch_data['incorrect_edits'],
            batch_data['scaffold_incorporated'],
            batch_data['prime_editor'],
            batch_data['replicate'],
            run_id,
            drug_id
        )
        data_points_added += 1
    
    print(f"\nRun '{run_name}' created successfully.")
    print(f"Added {data_points_added} data points to the database.")
    
    conn.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
