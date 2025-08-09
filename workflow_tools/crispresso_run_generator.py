#!/usr/bin/env python3
import os
import re
import sqlite3
import argparse
import glob
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def read_file_content(file_path: str) -> str:
    """Read and return the content of a file."""
    with open(file_path, 'r') as file:
        return file.read().strip()


def get_experiments_from_db(db_path: str) -> List[Dict]:
    """Retrieve all experiments from the database with relevant information."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    query = """
    SELECT id, name, variant, chromosome, genomic_location, edit, date
    FROM experiments
    ORDER BY date DESC, name
    """
    
    cursor.execute(query)
    experiments = [dict(row) for row in cursor.fetchall()]
    
    conn.close()
    return experiments


def display_experiments(experiments: List[Dict]) -> None:
    """Display experiments in a numbered list for user selection."""
    print("\nAvailable experiments:")
    print("-" * 80)
    print(f"{'#':<4} {'Name':<15} {'Variant':<10} {'Chr':<5} {'Location':<12} {'Edit':<20} {'Date':<12}")
    print("-" * 80)
    
    for i, exp in enumerate(experiments, 1):
        print(f"{i:<4} {exp['name']:<15} {exp['variant']:<10} {exp['chromosome']:<5} "
              f"{exp['genomic_location']:<12} {exp['edit']:<20} {exp['date']:<12}")


def get_experiment_entries(db_path: str, experiment_id: str) -> List[Dict]:
    """Get all experiment entries for a specific experiment."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    query = """
    SELECT id, name, pbs, rtt
    FROM experiment_entries
    WHERE experiment_id = ?
    """
    
    cursor.execute(query, (experiment_id,))
    entries = [dict(row) for row in cursor.fetchall()]
    
    conn.close()
    return entries


def get_protospacer_sequence(db_path: str, experiment_id: str) -> str:
    """Get the protospacer sequence for a specific experiment."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    query = """
    SELECT sense
    FROM protospacers
    WHERE experiment_id = ?
    """
    
    cursor.execute(query, (experiment_id,))
    result = cursor.fetchone()
    
    conn.close()
    
    if result:
        clean_extension =  re.sub(r'^[a-z]+|[a-z]+$', '', result[0])
        return clean_extension
    return ""


def get_extension_sequence(db_path: str, experiment_entry_id: str) -> str:
    """Get the extension sequence for a specific experiment entry."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    query = """
    SELECT sense
    FROM extensions
    WHERE experiment_entry_id = ?
    """
    
    cursor.execute(query, (experiment_entry_id,))
    result = cursor.fetchone()
    
    conn.close()
    
    if result:
        clean_extension = re.sub('^[a-z]+', '', result[0])
        return clean_extension
    return ""


def find_sequencing_files(seq_folder: str) -> Dict[str, Dict[str, str]]:
    """
    Find all sequencing files in the given folder and group them by sample name.
    Return a dictionary mapping sample names to their R1 and R2 files.
    """
    # Pattern to match sequencing files: {digits}_{digits}-{sample_name}_{letter+digits}_{R1 or R2}_001.fastq.gz
    pattern = re.compile(r'\d+_\d+-(.+?)_[A-Za-z0-9]+_(R[12])_001\.fastq\.gz')
    
    file_pairs = {}
    for file_path in glob.glob(os.path.join(seq_folder, "*.fastq.gz")):
        file_name = os.path.basename(file_path)
        match = pattern.match(file_name)
        
        if match:
            sample_name = match.group(1)
            read_type = match.group(2)
            
            if sample_name not in file_pairs:
                file_pairs[sample_name] = {}
                
            file_pairs[sample_name][read_type] = file_name
    
    # Only keep samples that have both R1 and R2 files
    return {sample: files for sample, files in file_pairs.items() 
            if 'R1' in files and 'R2' in files}


def validate_sample_names(sample_names: List[str]) -> Dict[str, Tuple[int, int]]:
    """
    Validate that sample names contain P{int} and R{int} and extract these values.
    Return a dictionary mapping sample names to (P value, R value) tuples.
    """
    p_pattern = re.compile(r'P(\d+)')
    r_pattern = re.compile(r'R(\d+)')
    
    valid_samples = {}
    for sample in sample_names:
        p_match = p_pattern.search(sample)
        r_match = r_pattern.search(sample)
        
        if p_match and r_match:
            p_value = int(p_match.group(1))
            r_value = int(r_match.group(1))
            valid_samples[sample] = (p_value, r_value)
    
    return valid_samples


def generate_output_file(output_path: str, 
                         sequencing_files: Dict[str, Dict[str, str]],
                         db_path: str, 
                         experiment_id: str,
                         amplicon_seq: str, 
                         scaffold_seq: str) -> None:
    """Generate the tab-separated output file with all the required information."""
    
    # Get protospacer sequence for the experiment
    protospacer_seq = get_protospacer_sequence(db_path, experiment_id)
    
    # Get experiment entries
    experiment_entries = get_experiment_entries(db_path, experiment_id)
    
    # Create a mapping from (PBS, RTT) to experiment entry ID
    entry_mapping = {}
    for entry in experiment_entries:
        entry_mapping[(entry['pbs'], entry['rtt'])] = entry['id']
    
    # Validate sample names and extract P and R values
    valid_samples = validate_sample_names(sequencing_files.keys())
    
    # Prepare output data
    output_data = []
    
    for sample_name, (p_value, r_value) in valid_samples.items():
        # Try to find an experiment entry with matching PBS and RTT values
        experiment_entry_id = entry_mapping.get((p_value, r_value))
        
        if experiment_entry_id:
            extension_seq = get_extension_sequence(db_path, experiment_entry_id)
            
            output_data.append({
                'name': sample_name,
                'fastq_r1': sequencing_files[sample_name]['R1'],
                'fastq_r2': sequencing_files[sample_name]['R2'],
                'prime_editing_pegRNA_extension_seq': extension_seq,
                'amplicon_seq': amplicon_seq,
                'prime_editing_pegRNA_scaffold_seq': scaffold_seq,
                'prime_editing_pegRNA_spacer_seq': protospacer_seq
            })
    
    # Write output to file
    with open(output_path, 'w') as out_file:
        # Write header
        header = [
            'name', 
            'fastq_r1', 
            'fastq_r2', 
            'prime_editing_pegRNA_extension_seq', 
            'amplicon_seq', 
            'prime_editing_pegRNA_scaffold_seq', 
            'prime_editing_pegRNA_spacer_seq'
        ]
        out_file.write('\t'.join(header) + '\n')
        
        # Write data
        for row in output_data:
            out_file.write('\t'.join([str(row[col]) for col in header]) + '\n')


def main():
    parser = argparse.ArgumentParser(description='Process sequencing data and generate a tab-separated output file.')
    parser.add_argument('--seq_folder', required=True, help='Folder containing sequencing data')
    parser.add_argument('--database', required=True, help='Path to the SQLite database')
    parser.add_argument('--amplicon_file', required=True, help='File containing the amplicon sequence')
    parser.add_argument('--scaffold_file', required=True, help='File containing the scaffold sequence')
    parser.add_argument('--output', default='output.txt', help='Output file path (default: output.txt)')
    
    args = parser.parse_args()
    
    # Read amplicon and scaffold sequences from files
    amplicon_seq = read_file_content(args.amplicon_file)
    scaffold_seq = read_file_content(args.scaffold_file)
    
    # Get experiments from database
    experiments = get_experiments_from_db(args.database)
    
    if not experiments:
        print("No experiments found in the database.")
        return
    
    # Display experiments and let user choose one
    display_experiments(experiments)
    
    while True:
        try:
            choice = int(input("\nSelect an experiment (1-{}): ".format(len(experiments))))
            if 1 <= choice <= len(experiments):
                selected_experiment = experiments[choice - 1]
                break
            else:
                print(f"Please enter a number between 1 and {len(experiments)}.")
        except ValueError:
            print("Please enter a valid number.")
    
    print(f"\nSelected experiment: {selected_experiment['name']} ({selected_experiment['variant']})")
    
    # Find sequencing files
    sequencing_files = find_sequencing_files(args.seq_folder)
    
    if not sequencing_files:
        print("No valid sequencing file pairs found in the specified folder.")
        return
    
    print(f"Found {len(sequencing_files)} sample(s) with valid R1/R2 file pairs.")
    
    # Generate output file
    generate_output_file(
        args.output,
        sequencing_files,
        args.database,
        selected_experiment['id'],
        amplicon_seq,
        scaffold_seq
    )
    
    print(f"\nOutput file successfully generated: {args.output}")


if __name__ == "__main__":
    main()
