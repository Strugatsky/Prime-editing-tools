#!/usr/bin/env python3
import sqlite3
import csv
import sys
import argparse

def convert_db_to_csv(db_file, output_file):
    """
    Convert SQLite database to CSV with specific formatting requirements.
    
    Args:
        db_file (str): Path to the SQLite database file
        output_file (str): Path to the output CSV file
    """
    # Connect to the database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Create CSV file
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write header
        header = ['experiment', 'run', 'prime_editor', 'PBS', 'RTT', 'replicate', 
                 'correct_edits', 'incorrect_edits', 'scaffold_incorporated']
        csv_writer.writerow(header)
        
        # SQL query to fetch all the required data
        query = """
        SELECT 
            e.name AS experiment,
            r.run_name AS run,
            d.prime_editor,
            ee.pbs AS PBS,
            ee.rtt AS RTT,
            d.replicate,
            d.correct_edits,
            d.incorrect_edits,
            d.scaffold_incorporated
        FROM data_points d
        JOIN experiment_entries ee ON d.experiment_entry_id = ee.id
        JOIN experiments e ON ee.experiment_id = e.id
        JOIN runs r ON d.run_id = r.id
        ORDER BY 
            e.name,
            r.run_name,
            d.prime_editor,
            ee.pbs,
            ee.rtt,
            d.replicate
        """
        
        # Execute query and write results to CSV
        cursor.execute(query)
        csv_writer.writerows(cursor.fetchall())
    
    # Close the connection
    conn.close()
    print(f"Data successfully exported to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Convert experiment database to CSV')
    parser.add_argument('db_file', help='SQLite database file')
    parser.add_argument('output_file', help='Output CSV file')
    
    args = parser.parse_args()
    
    try:
        convert_db_to_csv(args.db_file, args.output_file)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
