import os
import sys
import sqlite3
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import pandas as pd
import colorsys
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class ExperimentVisualizer:
    def __init__(self, root, db_path, output_dir, combined_plots=False):
        self.root = root
        self.root.title("Experiment Data Visualizer")
        self.root.geometry("800x600")
        
        # Initialize database and output directory
        self.db_path = db_path
        self.output_dir = output_dir
        self.combined_plots = combined_plots
        
        # Ensure output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Initialize variables
        self.experiments = []
        self.selected_experiment = None
        self.runs = []
        self.selected_run = None
        self.available_pbs_rtt = []
        self.selected_cells = set()
        self.dragging = False
        self.last_cell = None
        
        # Create GUI
        self.create_gui()
        
        # Load experiments
        self.load_experiments()
    
    def create_gui(self):
        # Create frames
        self.selection_frame = ttk.Frame(self.root, padding=10)
        self.selection_frame.pack(fill="both", expand=False)
        
        self.grid_frame = ttk.Frame(self.root, padding=10)
        self.grid_frame.pack(fill="both", expand=True)
        
        # Experiment selection
        ttk.Label(self.selection_frame, text="Select Experiment:").grid(row=0, column=0, sticky="w")
        self.experiment_var = tk.StringVar()
        self.experiment_dropdown = ttk.Combobox(self.selection_frame, textvariable=self.experiment_var, state="readonly", width=50)
        self.experiment_dropdown.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        self.experiment_dropdown.bind("<<ComboboxSelected>>", self.on_experiment_selected)
        
        # Run selection
        ttk.Label(self.selection_frame, text="Select Run:").grid(row=1, column=0, sticky="w")
        self.run_var = tk.StringVar()
        self.run_dropdown = ttk.Combobox(self.selection_frame, textvariable=self.run_var, state="readonly", width=50)
        self.run_dropdown.grid(row=1, column=1, padx=5, pady=5, sticky="ew")
        self.run_dropdown.bind("<<ComboboxSelected>>", self.on_run_selected)
        
        # PBS/RTT range labels
        ttk.Label(self.grid_frame, text="PBS/RTT Grid Selection:").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Label(self.grid_frame, text="PBS →").grid(row=1, column=2)
        ttk.Label(self.grid_frame, text="RTT\n↓").grid(row=2, column=0)
        
        # Plot button
        self.plot_button = ttk.Button(self.grid_frame, text="Generate Plots", command=self.generate_plots)
        self.plot_button.grid(row=3, column=1, pady=20)
        self.plot_button.state(['disabled'])
        
        # Clear selection button
        self.clear_button = ttk.Button(self.grid_frame, text="Clear Selection", command=self.clear_selection)
        self.clear_button.grid(row=3, column=2, pady=20)
        self.clear_button.state(['disabled'])
        
        # Status label
        self.status_var = tk.StringVar()
        self.status_var.set("Please select an experiment and run")
        self.status_label = ttk.Label(self.root, textvariable=self.status_var, foreground="blue")
        self.status_label.pack(pady=5)
    
    def load_experiments(self):
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Get all experiments with their names and IDs
            cursor.execute("SELECT id, name FROM experiments ORDER BY name")
            results = cursor.fetchall()
            
            # Format for dropdown: "1. Experiment Name"
            self.experiments = results
            experiment_options = [f"{i+1}. {name}" for i, (_, name) in enumerate(results)]
            
            self.experiment_dropdown['values'] = experiment_options
            
            conn.close()
            
            self.status_var.set(f"Loaded {len(results)} experiments")
        except Exception as e:
            self.status_var.set(f"Error loading experiments: {str(e)}")
    
    def on_experiment_selected(self, event):
        # Get the selected experiment index
        selection = self.experiment_dropdown.current()
        if selection >= 0:
            self.selected_experiment = self.experiments[selection][0]
            self.load_runs()
    
    def load_runs(self):
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Get runs for the selected experiment
            cursor.execute("""
                SELECT id, run_name 
                FROM runs 
                WHERE experiment_id = ? 
                ORDER BY run_name
            """, (self.selected_experiment,))
            
            results = cursor.fetchall()
            self.runs = results
            
            # Format for dropdown: "1. Run Name"
            run_options = [f"{i+1}. {name}" for i, (_, name) in enumerate(results)]
            
            self.run_dropdown['values'] = run_options
            self.run_var.set("")  # Clear selection
            self.status_var.set(f"Loaded {len(results)} runs for the selected experiment")
            
            conn.close()
        except Exception as e:
            self.status_var.set(f"Error loading runs: {str(e)}")
    
    def on_run_selected(self, event):
        # Get the selected run index
        selection = self.run_dropdown.current()
        if selection >= 0:
            self.selected_run = self.runs[selection][0]
            self.load_pbs_rtt_combinations()
    
    def load_pbs_rtt_combinations(self):
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Get valid PBS/RTT combinations for this run
            cursor.execute("""
                SELECT DISTINCT ee.pbs, ee.rtt
                FROM data_points dp
                JOIN experiment_entries ee ON dp.experiment_entry_id = ee.id
                WHERE dp.run_id = ?
                ORDER BY ee.pbs, ee.rtt
            """, (self.selected_run,))
            
            self.available_pbs_rtt = cursor.fetchall()
            
            if not self.available_pbs_rtt:
                self.status_var.set("No PBS/RTT combinations found for this run")
                return
            
            # Find min/max values for PBS and RTT
            pbs_values = [p for p, _ in self.available_pbs_rtt]
            rtt_values = [r for _, r in self.available_pbs_rtt]
            
            min_pbs, max_pbs = min(pbs_values), max(pbs_values)
            min_rtt, max_rtt = min(rtt_values), max(rtt_values)
            
            # Store current min/max values
            self.min_pbs = min_pbs
            self.max_pbs = max_pbs
            self.min_rtt = min_rtt
            self.max_rtt = max_rtt
            
            # Create grid
            self.create_grid()
            
            # Enable buttons
            self.plot_button.state(['!disabled'])
            self.clear_button.state(['!disabled'])
            
            self.status_var.set(f"Found PBS range: {min_pbs}-{max_pbs}, RTT range: {min_rtt}-{max_rtt}")
            
            conn.close()
        except Exception as e:
            self.status_var.set(f"Error loading PBS/RTT combinations: {str(e)}")
    
    def get_gradient_color(self, row, col, num_rows, num_cols):
        # Calculate position along diagonal (0 to 1)
        diagonal_pos = (row / num_rows + col / num_cols) / 2
        
        # Create a gradient through the HSV color space
        # Hue varies from 200 (blue) to 270 (purple)
        hue = 235 + diagonal_pos * 20
        saturation = 0.8
        value = 0.9 - diagonal_pos * 0.1
        
        # Convert HSV to RGB
        rgb = colorsys.hsv_to_rgb(hue/360, saturation, value)
        
        # Convert RGB to hex color string
        return f'#{int(rgb[0]*255):02x}{int(rgb[1]*255):02x}{int(rgb[2]*255):02x}'

    def create_grid(self):
        # Clear existing grid
        for widget in self.grid_frame.winfo_children():
            if isinstance(widget, tk.Canvas):
                widget.destroy()
        
        # Number of cells
        num_rows = self.max_rtt - self.min_rtt + 1
        num_cols = self.max_pbs - self.min_pbs + 1
        self.num_rows = num_rows
        self.num_cols = num_cols
        
        # Create canvas for grid
        canvas_size = 360  # Increased size for better visibility
        margin = 30  # Add margin for row/column numbers
        total_size = canvas_size + 2 * margin
        
        self.canvas = tk.Canvas(self.grid_frame, width=total_size, height=total_size)
        self.canvas.grid(row=2, column=1, padx=20, pady=20)
        
        # Calculate cell size
        self.cell_size = min(canvas_size / num_rows, canvas_size / num_cols)
        adjusted_width = self.cell_size * num_cols
        adjusted_height = self.cell_size * num_rows
        
        # Draw grid
        for i in range(num_rows + 1):
            # Horizontal lines
            y = margin + i * self.cell_size
            self.canvas.create_line(margin, y, margin + adjusted_width, y)
            
            if i < num_rows:
                rtt_value = self.max_rtt - i  # Reverse order so smaller RTT is at the bottom
                # Row numbers (left)
                self.canvas.create_text(margin/2, y + self.cell_size/2,
                                        text=str(rtt_value), anchor="e", font=("Arial", 8))
                # Row numbers (right)
                self.canvas.create_text(margin + adjusted_width + margin/2, y + self.cell_size/2,
                                        text=str(rtt_value), anchor="w", font=("Arial", 8))
        
        for i in range(num_cols + 1):
            # Vertical lines
            x = margin + i * self.cell_size
            self.canvas.create_line(x, margin, x, margin + adjusted_height)
            
            if i < num_cols:
                pbs_value = self.min_pbs + i
                # Column numbers (top)
                self.canvas.create_text(x + self.cell_size/2, margin/2,
                                        text=str(pbs_value), anchor="s", font=("Arial", 8))
                # Column numbers (bottom)
                self.canvas.create_text(x + self.cell_size/2, margin + adjusted_height + margin/2,
                                        text=str(pbs_value), anchor="n", font=("Arial", 8))
        
        # Shade cells that don't have data
        for row in range(num_rows):
            for col in range(num_cols):
                rtt_value = self.max_rtt - row  # Reverse order for RTT
                pbs_value = self.min_pbs + col
                
                if (pbs_value, rtt_value) not in self.available_pbs_rtt:
                    x1 = margin + col * self.cell_size
                    y1 = margin + row * self.cell_size
                    x2 = x1 + self.cell_size
                    y2 = y1 + self.cell_size
                    self.canvas.create_rectangle(x1, y1, x2, y2, fill="#f0f0f0")
        
        # Store margin for mouse event calculations
        self.margin = margin
        
        # Bind mouse events
        self.canvas.bind("<Button-1>", self.start_selection)
        self.canvas.bind("<B1-Motion>", self.update_selection)
        self.canvas.bind("<ButtonRelease-1>", self.end_selection)
        
        # Clear selection
        self.selected_cells.clear()

    def start_selection(self, event):
        self.dragging = True
        self.update_selection(event)

    def update_selection(self, event):
        if not self.dragging:
            return
            
        # Adjust for margin in cell coordinates calculation
        col = int((event.x - self.margin) // self.cell_size)
        row = int((event.y - self.margin) // self.cell_size)
        
        if 0 <= row < self.num_rows and 0 <= col < self.num_cols:
            # Convert grid coordinates to actual PBS/RTT values
            rtt_value = self.max_rtt - row  # Reverse order for RTT
            pbs_value = self.min_pbs + col
            
            # Only allow selection of valid PBS/RTT combinations
            if (pbs_value, rtt_value) in self.available_pbs_rtt:
                cell = (row, col)
                if cell != self.last_cell:
                    self.last_cell = cell
                    if cell in self.selected_cells:
                        self.selected_cells.remove(cell)
                    else:
                        self.selected_cells.add(cell)
                    self.draw_selection()

    def end_selection(self, event):
        self.dragging = False
        self.last_cell = None

    def draw_selection(self):
        self.canvas.delete("highlight")
        for row, col in self.selected_cells:
            x1 = self.margin + col * self.cell_size
            y1 = self.margin + row * self.cell_size
            x2 = x1 + self.cell_size
            y2 = y1 + self.cell_size
            color = self.get_gradient_color(row, col, self.num_rows, self.num_cols)
            self.canvas.create_rectangle(x1, y1, x2, y2, 
                                        fill=color, tags="highlight")
        self.canvas.tag_lower("highlight")

    def clear_selection(self):
        self.selected_cells.clear()
        self.canvas.delete("highlight")
        self.status_var.set("Selection cleared")

    def generate_plots(self):
        if not self.selected_cells:
            self.status_var.set("No cells selected. Please select at least one PBS/RTT combination.")
            return
    
        # Convert grid coordinates to actual PBS/RTT values
        selected_pbs_rtt = []
        for row, col in self.selected_cells:
            rtt_value = self.max_rtt - row  # Reverse order for RTT
            pbs_value = self.min_pbs + col
            selected_pbs_rtt.append((pbs_value, rtt_value))
    
        # Sort by PBS then RTT
        selected_pbs_rtt.sort()
    
        try:
            conn = sqlite3.connect(self.db_path)
        
            # Get data for selected PBS/RTT combinations
            # Using pandas for easier data manipulation
            query = """
                SELECT 
                    ee.pbs, 
                    ee.rtt, 
                    dp.correct_edits, 
                    dp.incorrect_edits, 
                    dp.scaffold_incorporated,
                    dp.prime_editor,
                    dp.replicate
                FROM data_points dp
                JOIN experiment_entries ee ON dp.experiment_entry_id = ee.id
                WHERE dp.run_id = ? AND (ee.pbs, ee.rtt) IN ({})
            """.format(','.join(['(?,?)'] * len(selected_pbs_rtt)))
        
            # Flatten the list of tuples for SQL parameters
            params = [self.selected_run] + [item for sublist in selected_pbs_rtt for item in sublist]
        
            df = pd.read_sql_query(query, conn, params=params)           
            if df.empty:
                self.status_var.set("No data found for the selected combinations")
                return
            
            # Create sample names
            df['sample'] = 'P' + df['pbs'].astype(str) + 'R' + df['rtt'].astype(str)
            
            # Get experiment and run name for plot titles
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM experiments WHERE id = ?", (self.selected_experiment,))
            experiment_name = cursor.fetchone()[0]
            
            cursor.execute("SELECT run_name FROM runs WHERE id = ?", (self.selected_run,))
            run_name = cursor.fetchone()[0]
            
            conn.close()
            
            # Generate plots
            self.plot_data(df, 'correct_edits', f"{experiment_name} - {run_name} - Correct Edits")
            self.plot_data(df, 'incorrect_edits', f"{experiment_name} - {run_name} - Incorrect Edits")
            self.plot_data(df, 'scaffold_incorporated', f"{experiment_name} - {run_name} - Scaffold Incorporated")
            
            self.status_var.set("Plots generated successfully")
            
        except Exception as e:
            self.status_var.set(f"Error generating plots: {str(e)}")
    
    def plot_data(self, df, metric, title):
        all_prime_editors = sorted(df['prime_editor'].unique())
        color_map = {}
        for editor in all_prime_editors:
            color_map[editor] = self.get_consistent_color(editor)

        if self.combined_plots:
            # Plot all prime editors on the same plot
            plt.figure(figsize=(10, 6), facecolor='white')
            plt.rcParams['axes.facecolor'] = 'white'

            for prime_editor in all_prime_editors:
                plt.figure(figsize=(10, 6), facecolor='white')
                plt.rcParams['axes.facecolor'] = 'white'
                
                for prime_editor in all_prime_editors:
                    if prime_editor not in df['prime_editor'].values:
                        continue  # Skip if this prime editor isn't in this dataset
                        
                    group_df = df[df['prime_editor'] == prime_editor]
                    
                    # Calculate average values per sample
                    avg_data = group_df.groupby('sample')[metric].mean().reset_index()
                    
                    # Sort by PBS first, then RTT
                    avg_data['pbs'] = avg_data['sample'].str.extract(r'P(\d+)R').astype(int)
                    avg_data['rtt'] = avg_data['sample'].str.extract(r'R(\d+)').astype(int)
                    avg_data = avg_data.sort_values(['pbs', 'rtt'])
                    
                    # Line plot of averages with consistent color for this prime editor
                    plt.plot(avg_data['sample'], avg_data[metric], marker='o', 
                             linewidth=2, markersize=8, color=color_map[prime_editor], 
                             label=prime_editor)
                    
                    # Plot individual replicates with matching color
                    for idx, sample in enumerate(avg_data['sample']):
                        sample_data = group_df[group_df['sample'] == sample]
                        x_jittered = [avg_data['sample'].tolist().index(sample) + np.random.normal(0, 0.05) 
                                     for _ in range(len(sample_data))]
                        plt.scatter(x_jittered, sample_data[metric], color=color_map[prime_editor], alpha=0.3)

            # Set x-axis to sample names from the first group (should be the same across all groups)
            all_samples = df.groupby('sample').first().reset_index()
            all_samples['pbs'] = all_samples['sample'].str.extract(r'P(\d+)R').astype(int)
            all_samples['rtt'] = all_samples['sample'].str.extract(r'R(\d+)').astype(int)
            all_samples = all_samples.sort_values(['pbs', 'rtt'])

            plt.xticks(range(len(all_samples)), all_samples['sample'], rotation=45)

            # Add labels, title and legend
            plt.ylabel(f"{metric.replace('_', ' ').title()} (%)")
            plt.xlabel("Sample")
            plt.title(f"{title} - All Prime Editors")
            plt.legend(title="Prime Editors")

            # Adjust layout
            plt.tight_layout()

            # Save plot with unique filename
            base_filename = f"{title.replace(' ', '_').replace('-', '_')}_All_Prime_Editors"
            filepath = os.path.join(self.output_dir, f"{base_filename}.png")

            # Check if file exists and create unique name if needed
            counter = 1
            while os.path.exists(filepath):
                filepath = os.path.join(self.output_dir, f"{base_filename}_{counter}.png")
                counter += 1

            plt.savefig(filepath, dpi=300, facecolor='white')
            plt.close()

        # Use a modern plot style with white background
        plt.figure(figsize=(10, 6), facecolor='white')
        plt.rcParams['axes.facecolor'] = 'white'
    
        # Group by prime editor and generate separate plots
        for prime_editor, group_df in df.groupby('prime_editor'):
            plt.figure(figsize=(10, 6), facecolor='white')
            plt.rcParams['axes.facecolor'] = 'white'
            
            # Calculate average values per sample
            avg_data = group_df.groupby('sample')[metric].mean().reset_index()
            
            # Sort by PBS first, then RTT
            avg_data['pbs'] = avg_data['sample'].str.extract(r'P(\d+)R').astype(int)
            avg_data['rtt'] = avg_data['sample'].str.extract(r'R(\d+)').astype(int)
            avg_data = avg_data.sort_values(['pbs', 'rtt'])
            
            # Line plot of averages with consistent color
            plt.plot(avg_data['sample'], avg_data[metric], marker='o', 
                     linewidth=2, markersize=8, color=color_map[prime_editor]) 
        
            # Plot individual replicates as grey dots
            for idx, sample in enumerate(avg_data['sample']):
                sample_data = group_df[group_df['sample'] == sample]
                x_jittered = [idx + np.random.normal(0, 0.05) for _ in range(len(sample_data))]
                plt.scatter(x_jittered, sample_data[metric], color='grey', alpha=0.5)
        
            # Set x-axis to show sample names
            plt.xticks(range(len(avg_data)), avg_data['sample'], rotation=45)
        
            # Add labels and title with prime editor info
            plt.ylabel(f"{metric.replace('_', ' ').title()} (%)")
            plt.xlabel("Sample")
            full_title = f"{title} - {prime_editor}"
            plt.title(full_title)
        
            # Adjust layout
            plt.tight_layout()
        
            # Save plot with unique filename to avoid overwriting
            base_filename = f"{full_title.replace(' ', '_').replace('-', '_')}"
            filepath = os.path.join(self.output_dir, f"{base_filename}.png")
        
            # Check if file exists and create unique name if needed
            counter = 1
            while os.path.exists(filepath):
                filepath = os.path.join(self.output_dir, f"{base_filename}_{counter}.png")
                counter += 1
        
            plt.savefig(filepath, dpi=300, facecolor='white')
            plt.close()

    def get_consistent_color(self, prime_editor):
        """Generate a consistent color based on the hash of the prime editor name"""
        import hashlib
        
        # Get hash value and convert to integer
        hash_value = int(hashlib.md5(prime_editor.encode()).hexdigest(), 16)
        
        # Use the hash to seed a random generator for consistent results
        import random
        rng = random.Random(hash_value)
        
        # Generate HSV color with fixed saturation and value, but varying hue
        h = rng.random()  # Hue between 0-1
        s = 0.7  # Saturation
        v = 0.9  # Value
        
        # Convert HSV to RGB
        import colorsys
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        
        return (r, g, b)


def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python script.py <database_path> <output_directory> [combined_plots]")
        sys.exit(1)
    
    db_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Optional argument for combined plots
    combined_plots = False
    if len(sys.argv) == 4 and sys.argv[3].lower() in ['true', 'yes', '1', 'combined']:
        combined_plots = True
    
    # Check if database exists
    if not os.path.isfile(db_path):
        print(f"Error: Database file {db_path} not found")
        sys.exit(1)
    
    # Create tkinter window
    root = tk.Tk()
    app = ExperimentVisualizer(root, db_path, output_dir, combined_plots)
    root.mainloop()

if __name__ == "__main__":
    main()
