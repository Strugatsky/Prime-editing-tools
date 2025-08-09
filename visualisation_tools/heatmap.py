import argparse
import sqlite3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


def connect_to_database(db_path):
    """Connect to the SQLite database."""
    if not os.path.exists(db_path):
        print(f"Error: Database file '{db_path}' not found.")
        exit(1)
    
    try:
        conn = sqlite3.connect(db_path)
        return conn
    except sqlite3.Error as e:
        print(f"Error connecting to database: {e}")
        exit(1)


def get_experiments(conn):
    """Get list of experiments from the database."""
    cursor = conn.cursor()
    query = """
        SELECT id, rowid, name, variant, chromosome, genomic_location, edit, date
        FROM experiments
        ORDER BY date DESC
    """
    
    try:
        cursor.execute(query)
        experiments = cursor.fetchall()
        return experiments
    except sqlite3.Error as e:
        print(f"Error fetching experiments: {e}")
        exit(1)


def display_experiments(experiments):
    """Display numbered list of experiments for user selection."""
    print("\nAvailable Experiments:")
    print("-" * 100)
    print(f"{'#':<5} {'Name':<20} {'Variant':<15} {'Chr':<5} {'Location':<15} {'Edit':<20} {'Date':<12}")
    print("-" * 100)

    for i, experiment in enumerate(experiments, 1):
        exp_id, rowid, name, variant, chromosome, location, edit, date = experiment
        print(f"{i:<5} {name:<20} {variant:<15} {chromosome:<5} {location:<15} {edit:<20} {date:<12}")


def get_experiment_data(conn, experiment_id, selected_runs=None, normalize=False):
    """Get all data for a specific experiment."""
    cursor = conn.cursor()
    
    # Get all prime editors for this experiment
    if selected_runs:
        query = """
        SELECT DISTINCT d.prime_editor
        FROM data_points d
        JOIN experiment_entries e ON d.experiment_entry_id = e.id
        WHERE e.experiment_id = ? AND d.run_id IN ({})
        """.format(','.join('?' * len(selected_runs)))
        cursor.execute(query, (experiment_id, *selected_runs))
    else:
        query = """
        SELECT DISTINCT d.prime_editor
        FROM data_points d
        JOIN experiment_entries e ON d.experiment_entry_id = e.id
        WHERE e.experiment_id = ?
        """
        cursor.execute(query, (experiment_id,))
    
    prime_editors = [row[0] for row in cursor.fetchall()]

    # Get all drugs used in this experiment
    if selected_runs:
        query = """
        SELECT DISTINCT dr.id, dr.name
        FROM data_points d
        JOIN experiment_entries e ON d.experiment_entry_id = e.id
        JOIN drugs dr ON d.drug_id = dr.id
        WHERE e.experiment_id = ? AND d.run_id IN ({})
        """.format(','.join('?' * len(selected_runs)))
        cursor.execute(query, (experiment_id, *selected_runs))
    else:
        query = """
        SELECT DISTINCT dr.id, dr.name
        FROM data_points d
        JOIN experiment_entries e ON d.experiment_entry_id = e.id
        JOIN drugs dr ON d.drug_id = dr.id
        WHERE e.experiment_id = ?
        """
        cursor.execute(query, (experiment_id,))

    drugs = {row[0]: row[1] for row in cursor.fetchall()}

    # If no drugs found (for compatibility with old databases), add a default "None" drug
    if not drugs:
        drugs = {"00000000-0000-0000-0000-000000000000": "None"}
    
    # Get all run IDs for this experiment
    query = """
    SELECT id, run_name
    FROM runs
    WHERE experiment_id = ?
    """
    cursor.execute(query, (experiment_id,))
    all_runs = {row[0]: row[1] for row in cursor.fetchall()}
    if selected_runs:
        runs = {run_id: name for run_id, name in all_runs.items() if run_id in selected_runs}
    else:
        runs = all_runs
    
    # Get experiment entries (PBS and RTT combinations)
    query = """
    SELECT id, pbs, rtt
    FROM experiment_entries
    WHERE experiment_id = ?
    """
    cursor.execute(query, (experiment_id,))
    entries = {row[0]: (row[1], row[2]) for row in cursor.fetchall()}
    
    # Get all data points
    if selected_runs:
        query = """
        SELECT d.experiment_entry_id, d.correct_edits, d.incorrect_edits, 
               d.scaffold_incorporated, d.prime_editor, d.replicate, d.run_id, 
               COALESCE(d.drug_id, '00000000-0000-0000-0000-000000000000') as drug_id
        FROM data_points d
        JOIN experiment_entries e ON d.experiment_entry_id = e.id
        WHERE e.experiment_id = ? AND d.run_id IN ({})
        """.format(','.join('?' * len(selected_runs)))
        cursor.execute(query, (experiment_id, *selected_runs))
    else:
        query = """
        SELECT d.experiment_entry_id, d.correct_edits, d.incorrect_edits, 
               d.scaffold_incorporated, d.prime_editor, d.replicate, d.run_id,
               COALESCE(d.drug_id, '00000000-0000-0000-0000-000000000000') as drug_id
        FROM data_points d
        JOIN experiment_entries e ON d.experiment_entry_id = e.id
        WHERE e.experiment_id = ?
        """
        cursor.execute(query, (experiment_id,))

    data_points = cursor.fetchall()
    
    # Collect data for heatmaps
    heatmap_data = defaultdict(lambda: defaultdict(list))
    
    # Find the common PBS/RTT combination for normalization if requested
    if normalize:
        # Track which PBS/RTT combinations are in which runs
        pbs_rtt_in_runs = defaultdict(set)
        for entry_id, _, _, _, _, _, run_id in data_points:
            if entry_id in entries:
                pbs, rtt = entries[entry_id]
                pbs_rtt_in_runs[(pbs, rtt)].add(run_id)
        
        # Find combinations present in all runs
        all_runs = set(runs.keys())
        common_pbs_rtt = [combo for combo, run_set in pbs_rtt_in_runs.items() 
                         if run_set == all_runs]
        
        if not common_pbs_rtt:
            print("Warning: No common PBS/RTT combination found across all runs. Normalization disabled.")
            normalize = False
        else:
            # Use the first common combination for normalization
            norm_pbs, norm_rtt = common_pbs_rtt[0]
            print(f"Normalizing data using PBS={norm_pbs}, RTT={norm_rtt} as reference point.")
    
    # Process data points
    run_normalization_factors = {}
    if normalize:
        # Calculate reference values for each run
        for run_id in runs:
            ref_values = []
            for entry_id, correct, incorrect, scaffold, editor, replicate, r_id in data_points:
                if r_id == run_id and entry_id in entries:
                    pbs, rtt = entries[entry_id]
                    if pbs == norm_pbs and rtt == norm_rtt:
                        # Store correct_edits for normalization (could use other metrics)
                        ref_values.append(correct)
            
            if ref_values:
                run_normalization_factors[run_id] = np.mean(ref_values)
    
    # Determine global normalization factor (average of reference points)
    if normalize and run_normalization_factors:
        global_factor = np.mean(list(run_normalization_factors.values()))
    
    # Process all data points
    for entry_id, correct, incorrect, scaffold, editor, replicate, run_id, drug_id in data_points:
        if entry_id in entries:
            pbs, rtt = entries[entry_id]
            
            # Apply normalization if enabled
            norm_factor = 1.0
            if normalize and run_id in run_normalization_factors and run_normalization_factors[run_id] > 0:
                norm_factor = global_factor / run_normalization_factors[run_id]
            
            # Store data for each metric and prime editor
            drug_name = drugs.get(drug_id, "Unknown")
            heatmap_data[(editor, drug_name, 'correct_edits')][(pbs, rtt)].append(correct * norm_factor)
            heatmap_data[(editor, drug_name, 'incorrect_edits')][(pbs, rtt)].append(incorrect * norm_factor)
            heatmap_data[(editor, drug_name, 'scaffold_incorporated')][(pbs, rtt)].append(scaffold * norm_factor)
    
    # Calculate averages for replicates
    averaged_data = {}
    for key, values in heatmap_data.items():
        averaged_data[key] = {pos: np.mean(vals) for pos, vals in values.items()}
    
    # Get experiment details for titles
    query = "SELECT name, variant FROM experiments WHERE id = ?"
    cursor.execute(query, (experiment_id,))
    exp_details = cursor.fetchone()
    exp_name = exp_details[0]
    exp_variant = exp_details[1] if exp_details[1] else ""
    
    return averaged_data, prime_editors, drugs, exp_name, exp_variant, runs


def create_heatmaps(data, prime_editors, drugs, exp_name, exp_variant, normalize, run_info=""):
    """Create and save heatmaps for each prime editor and metric."""
    metrics = ['correct_edits', 'incorrect_edits', 'scaffold_incorporated']
    
    # Create directory for output if it doesn't exist
    output_dir = f"heatmaps_{exp_name.replace(' ', '_')}"
    if exp_variant:
        output_dir += f"_{exp_variant.replace(' ', '_')}"
    if run_info:
        output_dir += f"_{run_info}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    norm_text = "_normalized" if normalize else ""
    
    for editor in prime_editors:
        for drug_name in drugs.values():
            for metric in metrics:
                key = (editor, drug_name, metric)
                if key not in data or not data[key]:
                    print(key, "key not found in data, continuing")
                    continue
                
                # Extract PBS and RTT values
                positions = list(data[key].keys())
                pbs_values = sorted(set([pos[0] for pos in positions]))
                rtt_values = sorted(set([pos[1] for pos in positions]))
                
                # Create matrix for heatmap
                matrix = np.full((len(pbs_values), len(rtt_values)), np.nan)
                
                # Fill matrix with data
                for i, pbs in enumerate(pbs_values):
                    for j, rtt in enumerate(rtt_values):
                        if (pbs, rtt) in data[key]:
                            matrix[i, j] = data[key][(pbs, rtt)]
                
                # Create heatmap
                plt.figure(figsize=(10, 8))
                mask = np.isnan(matrix)
                ax = sns.heatmap(matrix, annot=True, fmt=".2f", cmap="viridis", 
                            xticklabels=rtt_values, yticklabels=pbs_values,
                                 mask=mask, annot_kws={"size": 9})
                
                metric_name = metric.replace('_', ' ').title()
                drug_display = f"Drug: {drug_name}" if drug_name != "None" else "No Drug"
                plt.title(f"{exp_name} - {exp_variant}\n{editor} - {drug_display} - {metric_name}")
                plt.xlabel("RTT Length")
                plt.ylabel("PBS Length")
                
                # Save as SVG
                filename = f"{output_dir}/{editor}_{metric}_{drug_name}{norm_text}.svg"
                plt.savefig(filename, format='svg')
                plt.close()
                
                print(f"Saved heatmap to {filename}")

def select_runs(conn, experiment_id, all_runs):
    """Display all runs for an experiment and let user select one or more."""
    runs = [(run_id, run_name) for run_id, run_name in all_runs.items()]

    if not runs:
        print("No runs found for this experiment.")
        return []

    print("\nAvailable Runs:")
    for i, (run_id, run_name) in enumerate(runs, 1):
        print(f"{i}. {run_name}")

    print("\nEnter run numbers to include (comma-separated, or 'all' for all runs):")
    selection = input("> ").strip()

    if selection.lower() == 'all':
        return [run_id for run_id, _ in runs]

    try:
        selected_indices = [int(idx.strip()) for idx in selection.split(',')]
        selected_runs = [runs[idx-1][0] for idx in selected_indices if 1 <= idx <= len(runs)]

        if not selected_runs:
            print("No valid runs selected. Using all runs.")
            return [run_id for run_id, _ in runs]

        return selected_runs
    except ValueError:
        print("Invalid selection. Using all runs.")
        return [run_id for run_id, _ in runs]

def get_all_runs(conn, experiment_id):
    """Get all runs for an experiment."""
    cursor = conn.cursor()
    query = """
    SELECT id, run_name
    FROM runs
    WHERE experiment_id = ?
    """
    cursor.execute(query, (experiment_id,))
    runs = {row[0]: row[1] for row in cursor.fetchall()}
    return runs

def main():
    parser = argparse.ArgumentParser(description='Generate heatmaps from experiment data.')
    parser.add_argument('database', help='Path to the SQLite database file')
    parser.add_argument('--normalize', action='store_true', help='Normalize data across runs')
    args = parser.parse_args()
    
    conn = connect_to_database(args.database)
    experiments = get_experiments(conn)
    
    if not experiments:
        print("No experiments found in the database.")
        conn.close()
        exit(0)
    
    display_experiments(experiments)
    
    try:
        choice = int(input("\nSelect an experiment (1-{}): ".format(len(experiments))))
        if choice < 1 or choice > len(experiments):
            print("Invalid choice.")
            conn.close()
            exit(1)
    except ValueError:
        print("Please enter a valid number.")
        conn.close()
        exit(1)
    
    experiment_id = experiments[choice-1][0]
    all_runs = get_all_runs(conn, experiment_id)
    selected_runs = select_runs(conn, experiment_id, all_runs)
    data, prime_editors, drugs, exp_name, exp_variant, runs = get_experiment_data(conn, experiment_id, selected_runs, args.normalize)    
    if selected_runs and len(selected_runs) < len(all_runs):
        selected_run_names = [all_runs[run_id] for run_id in selected_runs]
        if len(selected_run_names) <= 2:
            run_info = "_" + "_".join(selected_run_names)
        else:
            run_info = f"_{len(selected_run_names)}_selected_runs"
    else:
        run_info = "" 
    create_heatmaps(data, prime_editors, drugs, exp_name, exp_variant, args.normalize, run_info)
    
    conn.close()
    print("\nAll heatmaps generated successfully.")


if __name__ == "__main__":
    main()
