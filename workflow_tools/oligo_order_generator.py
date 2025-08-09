#!/usr/bin/env python3
import csv
import sys
import argparse
import sqlite3
import os
import re

def validate_and_replace_scaffold(sequence):
    if sequence.startswith('gtgc'):
        return 'gtcc' + sequence[4:]
    raise ValueError("sequence doesn't start with gtgc")

def select_experiment(db_path):
    """
    Connect to the database and let the user select an experiment.
    
    Args:
        db_path (str): Path to SQLite database
    
    Returns:
        tuple: (experiment_id, experiment_name, experiment_variant)
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Fetch all experiments
        cursor.execute("""
            SELECT id, name, variant, chromosome, genomic_location, date
            FROM experiments
            ORDER BY date DESC
        """)
        
        experiments = cursor.fetchall()
        
        if not experiments:
            print("No experiments found in the database.")
            conn.close()
            sys.exit(1)
        
        print("\nAvailable experiments:")
        print("----------------------")
        for i, exp in enumerate(experiments, 1):
            exp_id, name, variant, chromosome, location, date = exp
            print(f"{i}. {name} (Variant: {variant}, Chr: {chromosome}, Loc: {location}, Date: {date}, ID: {exp_id})")
        
        while True:
            try:
                choice = int(input("\nSelect an experiment (1-{}): ".format(len(experiments))))
                if 1 <= choice <= len(experiments):
                    selected = experiments[choice-1]
                    print(f"\nSelected: {selected[1]} (Variant: {selected[2]})")
                    return selected[0], selected[1], selected[2]
                else:
                    print("Invalid choice. Please try again.")
            except ValueError:
                print("Please enter a number.")
    finally:
        conn.close()

def setup_database_tables(db_path):
    """
    Set up new tables in the database if they don't exist.
    
    Args:
        db_path (str): Path to SQLite database
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Check if protospacers table exists, create if not
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS protospacers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                experiment_id TEXT NOT NULL,
                sense TEXT NOT NULL,
                antisense TEXT NOT NULL,
                FOREIGN KEY (experiment_id) REFERENCES experiments(id)
            )
        """)
        
        # Check if extensions table exists, create if not
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS extensions (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                experiment_entry_id INTEGER NOT NULL,
                sense TEXT NOT NULL,
                antisense TEXT NOT NULL,
                FOREIGN KEY (experiment_entry_id) REFERENCES experiment_entries(id)
            )
        """)
        
        # Add columns to experiments table if they don't exist
        cursor.execute("PRAGMA table_info(experiments)")
        columns = [info[1] for info in cursor.fetchall()]
        
        if 'edit_position' not in columns:
            cursor.execute("ALTER TABLE experiments ADD COLUMN edit_position TEXT")
        if 'pam' not in columns:
            cursor.execute("ALTER TABLE experiments ADD COLUMN pam TEXT")
        if 'pam_strand' not in columns:
            cursor.execute("ALTER TABLE experiments ADD COLUMN pam_strand TEXT")
            
        # Add name column to experiment_entries if it doesn't exist
        cursor.execute("PRAGMA table_info(experiment_entries)")
        columns = [info[1] for info in cursor.fetchall()]
        
        if 'name' not in columns:
            cursor.execute("ALTER TABLE experiment_entries ADD COLUMN name TEXT")
        if 'score' not in columns:
            cursor.execute("ALTER TABLE experiment_entries ADD COLUMN score TEXT")
            
        conn.commit()
    finally:
        conn.close()

def process_csv_and_update_db(input_file, output_file, db_path, prefix=None, use_scaffold_2=False):
    """
    Process input CSV, update database, and create output CSV.
    
    Args:
        input_file (str): Path to input CSV file
        output_file (str): Path to output CSV file
        db_path (str): Path to SQLite database
        prefix (str, optional): Prefix to add to variant names
        use_scaffold_2 (bool): Whether to replace 'gtgc' with 'gtcc' in sense sequences
    """
    # First, let the user select an experiment
    experiment_id, experiment_name, variant = select_experiment(db_path)
    
    # Set up database tables
    setup_database_tables(db_path)
    
    # Connect to the database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Process the CSV file
        with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            # Read input CSV with headers
            reader = csv.DictReader(infile)
            writer = csv.writer(outfile)
            
            # Get the first row to extract common data
            rows = list(reader)
            if not rows:
                print("Input CSV file is empty.")
                return
            
            first_row = rows[0]
            
            # Update experiment table with common data
            cursor.execute("""
                UPDATE experiments
                SET edit_position = ?, pam = ?, pam_strand = ?
                WHERE id = ?
            """, (first_row['EditPos.'], first_row['PAM'], first_row['PAM.Strand'], experiment_id))
            
            # Add protospacer data (should be the same for each experiment)
            cursor.execute("""
                INSERT OR REPLACE INTO protospacers (experiment_id, sense, antisense)
                VALUES (?, ?, ?)
            """, (experiment_id, first_row['Protospacer.Sense.'], first_row['Protospacer.Antisense.']))
            
            # Process each row
            for row_number, row in enumerate(rows, start=1):
                try:
                    # Extract needed values
                    pbs = row['PBS']
                    rtt = row['RTT']
                    score = row['Score']
                    ext_sense = row['Extension.Sense.']
                    ext_antisense = row['Extension.Antisense.']
                    
                    # Process sense sequence if scaffold 2 is requested
                    if use_scaffold_2:
                        try:
                            ext_sense = validate_and_replace_scaffold(ext_sense)
                        except ValueError as e:
                            raise ValueError(f"Row {row_number} - {str(e)}")
                    
                    # Create base identifier
                    base_name = f"{prefix + '_' if prefix else ''}{variant}_P{pbs}_R{rtt}"
                    
                    # Find the corresponding experiment_entry
                    cursor.execute("""
                        SELECT id FROM experiment_entries
                        WHERE experiment_id = ? AND pbs = ? AND rtt = ?
                    """, (experiment_id, pbs, rtt))
                    
                    entry = cursor.fetchone()
                    if entry:
                        entry_id = entry[0]
                        
                        # Update the experiment_entry with name and score
                        cursor.execute("""
                            UPDATE experiment_entries
                            SET name = ?, score = ?
                            WHERE id = ?
                        """, (base_name, score, entry_id))
                        
                        # Add extension data
                        cursor.execute("""
                            INSERT INTO extensions (experiment_entry_id, sense, antisense)
                            VALUES (?, ?, ?)
                        """, (entry_id, ext_sense, ext_antisense))
                        
                        # Write to output CSV
                        writer.writerow([f"{base_name}_S", ext_sense])
                        writer.writerow([f"{base_name}_AS", ext_antisense])
                    else:
                        print(f"Warning: No matching experiment entry found for PBS={pbs}, RTT={rtt}")
                
                except ValueError as e:
                    print(f"Error processing row {row_number}: {str(e)}", file=sys.stderr)
                    conn.rollback()
                    sys.exit(1)
            
            # Commit changes
            conn.commit()
            print(f"Database updated successfully for experiment: {experiment_name}")
    finally:
        conn.close()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process CSV file and update SQLite database for variant analysis')
    parser.add_argument('input_file', help='Input CSV file path')
    parser.add_argument('output_file', help='Output CSV file path')
    parser.add_argument('db_file', help='SQLite database file path')
    parser.add_argument('--prefix', help='Optional prefix for variant names', default=None)
    parser.add_argument('--use_scaffold_2', action='store_true', 
                        help='Replace gtgc with gtcc in sense sequences')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if files exist
    if not os.path.isfile(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isfile(args.db_file):
        print(f"Error: Database file '{args.db_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    try:
        process_csv_and_update_db(
            args.input_file, 
            args.output_file, 
            args.db_file, 
            args.prefix, 
            args.use_scaffold_2
        )
        print(f"Processing complete. Output written to {args.output_file}")
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
