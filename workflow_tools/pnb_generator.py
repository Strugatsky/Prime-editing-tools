import tkinter as tk
import ttkbootstrap as ttk
import csv
from datetime import datetime
import colorsys
from tkinter import messagebox, filedialog
import os
import sqlite3
import uuid

class GridSelectorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Grid Selector")
        self.root.geometry("1200x645")  # Wider window, shorter height

        # Variables
        self.selected_cells = set()
        self.dragging = False
        self.last_cell = None
        self.num_cells = 0  # Initialize num_cells
        self.output_path = tk.StringVar(value="pnb_batch.csv")
        self.db_path = tk.StringVar(value="experiments.db")
        self.experiment_name = tk.StringVar(value="Experiment_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
        self.db_connection = None
        
        # Create main container frames
        self.left_frame = ttk.Frame(self.root, padding="10")
        self.left_frame.pack(side="left", fill="both", expand=False, padx=10, pady=10)
        
        self.right_frame = ttk.Frame(self.root, padding="10")
        self.right_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)
        
        # Create frames in the left column
        self.create_range_frame()
        self.create_input_frame()
        self.create_output_file_frame()
        self.create_database_frame()
        
        # Create frames in the right column
        self.create_grid_frame()
        self.create_button_frame()
        
        # Initial grid creation
        self.create_grid()

    def create_range_frame(self):
        range_frame = ttk.LabelFrame(self.left_frame, text="Grid Range", padding="10")
        range_frame.pack(fill="x", pady=5)
        
        ttk.Label(range_frame, text="Min Value:").grid(row=0, column=0, padx=5)
        self.min_var = tk.StringVar(value="7")
        ttk.Entry(range_frame, textvariable=self.min_var, width=10).grid(row=0, column=1, padx=5)
        
        ttk.Label(range_frame, text="Max Value:").grid(row=0, column=2, padx=5)
        self.max_var = tk.StringVar(value="16")
        ttk.Entry(range_frame, textvariable=self.max_var, width=10).grid(row=0, column=3, padx=5)
        
        ttk.Button(range_frame, text="Update Grid", command=self.create_grid, bootstyle="info").grid(row=0, column=4, padx=5)

    def create_grid_frame(self):
        self.grid_frame = ttk.LabelFrame(self.right_frame, text="Selection Grid", padding="10")
        self.grid_frame.pack(fill="both", expand=True, pady=5)
        
        # Labels for axes
        ttk.Label(self.grid_frame, text="RTT length (nts)", anchor="center").grid(row=0, column=0, columnspan=2)
        ttk.Label(self.grid_frame, text="PBS length (nts)", anchor="center").grid(row=1, column=0, rowspan=2)

    def create_input_frame(self):
        input_frame = ttk.LabelFrame(self.left_frame, text="Information", padding="10")
        input_frame.pack(fill="x", pady=5)
        
        self.inputs = {}
        fields = ["Variant", "Chromosome", "Genomic Location", "Edit", "Gene Orientation"]
        defaults = ["HEK3_1CTTins", "9", "107422356", "insCTT", "+"]
        
        for i, field in enumerate(fields):
            ttk.Label(input_frame, text=f"{field}:").grid(row=i, column=0, padx=5, pady=2, sticky="e")
            self.inputs[field] = ttk.Entry(input_frame, width=30)
            self.inputs[field].insert(tk.END, defaults[i])
            self.inputs[field].grid(row=i, column=1, padx=5, pady=2, sticky="w")

    def create_output_file_frame(self):
        file_frame = ttk.LabelFrame(self.left_frame, text="Output File", padding="10")
        file_frame.pack(fill="x", pady=5)

        # File path entry
        ttk.Label(file_frame, text="Save as:").grid(row=0, column=0, padx=5)
        self.file_entry = ttk.Entry(file_frame, textvariable=self.output_path, width=30)
        self.file_entry.grid(row=0, column=1, padx=5, sticky="ew")

        # Browse button
        ttk.Button(file_frame, text="Browse", command=self.browse_output_file).grid(row=0, column=2, padx=5)

        # Configure grid weights
        file_frame.columnconfigure(1, weight=1)

    def create_database_frame(self):
        db_frame = ttk.LabelFrame(self.left_frame, text="Database Configuration", padding="10")
        db_frame.pack(fill="x", pady=5)

        # Experiment name
        ttk.Label(db_frame, text="Experiment Name:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        ttk.Entry(db_frame, textvariable=self.experiment_name, width=30).grid(row=0, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
        
        # Database path
        ttk.Label(db_frame, text="Database File:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        ttk.Entry(db_frame, textvariable=self.db_path, width=30).grid(row=1, column=1, padx=5, pady=5, sticky="ew")
        ttk.Button(db_frame, text="Browse", command=self.browse_database_file).grid(row=1, column=2, padx=5, pady=5)

        # Database actions
        ttk.Button(db_frame, text="Connect", command=self.connect_database, bootstyle="primary").grid(row=2, column=0, padx=5, pady=5)
        ttk.Button(db_frame, text="Save to DB", command=self.save_to_database, bootstyle="success").grid(row=2, column=1, padx=5, pady=5, sticky="ew")
        
        # Configure grid weights
        db_frame.columnconfigure(1, weight=1)

    def browse_output_file(self):
        # Get the current directory and filename
        current_dir = os.path.dirname(self.output_path.get()) or os.getcwd()
        current_file = os.path.basename(self.output_path.get())
        
        # Open file dialog
        filename = filedialog.asksaveasfilename(
            initialdir=current_dir,
            initialfile=current_file,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            self.output_path.set(filename)

    def browse_database_file(self):
        # Get the current directory and filename
        current_dir = os.path.dirname(self.db_path.get()) or os.getcwd()
        current_file = os.path.basename(self.db_path.get())
        
        # Open file dialog
        filename = filedialog.asksaveasfilename(
            initialdir=current_dir,
            initialfile=current_file,
            defaultextension=".db",
            filetypes=[("SQLite Database", "*.db"), ("All files", "*.*")]
        )
        
        if filename:
            self.db_path.set(filename)

    def create_button_frame(self):
        button_frame = ttk.Frame(self.right_frame, padding="10")
        button_frame.pack(fill="x", pady=5)
        
        ttk.Button(button_frame, text="Generate CSV", bootstyle="success", command=self.generate_csv).pack(side="right", padx=5)
        ttk.Button(button_frame, text="Clear Selection", bootstyle="danger", command=self.clear_selection).pack(side="right", padx=5)

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
        
        # Get range values
        try:
            min_val = int(self.min_var.get())
            max_val = int(self.max_var.get())
            self.num_cells = max_val - min_val + 1  # Store num_cells as instance variable
        except ValueError:
            return
        
        # Create canvas for grid
        canvas_size = 360  # Increased size for better visibility
        margin = 30  # Add margin for row/column numbers
        total_size = canvas_size + 2 * margin
        
        self.canvas = tk.Canvas(self.grid_frame, width=total_size, height=total_size)
        self.canvas.grid(row=2, column=1, padx=20, pady=20)
        
        # Calculate cell size
        self.cell_size = canvas_size / self.num_cells
        
        # Draw grid
        for i in range(self.num_cells + 1):
            # Draw lines
            x = margin + i * self.cell_size
            y = margin + i * self.cell_size
            
            # Vertical and horizontal grid lines
            self.canvas.create_line(x, margin, x, margin + canvas_size)
            self.canvas.create_line(margin, y, margin + canvas_size, y)
            
            if i < self.num_cells:
                value = min_val + i
                # Column numbers (top)
                self.canvas.create_text(x + self.cell_size/2, margin/2,
                                      text=str(value), anchor="s", font=("Arial", 8))
                # Column numbers (bottom)
                self.canvas.create_text(x + self.cell_size/2, margin + canvas_size + margin/2,
                                      text=str(value), anchor="n", font=("Arial", 8))
                # Row numbers (left)
                self.canvas.create_text(margin/2, y + self.cell_size/2,
                                      text=str(value), anchor="e", font=("Arial", 8))
                # Row numbers (right)
                self.canvas.create_text(margin + canvas_size + margin/2, y + self.cell_size/2,
                                      text=str(value), anchor="w", font=("Arial", 8))
        
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
        
        if 0 <= row < self.num_cells and 0 <= col < self.num_cells:
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
            color = self.get_gradient_color(row, col, self.num_cells, self.num_cells)
            self.canvas.create_rectangle(x1, y1, x2, y2, 
                                      fill=color, tags="highlight")
        self.canvas.tag_lower("highlight")

    def clear_selection(self):
        self.selected_cells.clear()
        self.canvas.delete("highlight")

    def generate_csv(self):
        if not self.selected_cells:
            messagebox.showwarning("Warning", "No cells selected")
            return
            
        # Prepare data
        min_val = int(self.min_var.get())
        headers = ["Variant", "Chromosome", "GenomicLocation", "Edit", 
                  "GeneOrientation", "PBS", "RTT"]        
        rows = []
        for row, col in self.selected_cells:
            row_data = [
                self.inputs["Variant"].get(),
                self.inputs["Chromosome"].get(),
                self.inputs["Genomic Location"].get(),
                self.inputs["Edit"].get(),
                self.inputs["Gene Orientation"].get(),
                str(min_val + col),  # PBS length
                str(min_val + row)   # RTT length
            ]
            rows.append(row_data)        
        # Write CSV
        try:
            with open(self.output_path.get(), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(rows)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")
            return
            
        messagebox.showinfo("Success", f"CSV file '{self.output_path.get()}' has been generated")

    # Database methods
    def connect_database(self):
        db_path = self.db_path.get()
        if not db_path:
            messagebox.showwarning("Warning", "Please specify a database file path")
            return
            
        try:
            # Close existing connection if any
            if self.db_connection:
                self.db_connection.close()
                
            # Create directory if it doesn't exist
            db_dir = os.path.dirname(db_path)
            if db_dir and not os.path.exists(db_dir):
                os.makedirs(db_dir)
                
            # Connect to database
            self.db_connection = sqlite3.connect(db_path)
            self.initialize_database()
            messagebox.showinfo("Success", f"Connected to database: {db_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to connect to database: {str(e)}")
            self.db_connection = None

    def initialize_database(self):
        """Initialize database tables if they don't exist"""
        cursor = self.db_connection.cursor()
        
        # Create experiments table
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS experiments (
            id TEXT PRIMARY KEY,
            name TEXT NOT NULL,
            date TEXT NOT NULL,
            variant TEXT,
            chromosome TEXT,
            genomic_location TEXT,
            edit TEXT,
            gene_orientation TEXT
        )
        ''')
        
        # Create entries table for PBS/RTT combinations
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS experiment_entries (
            id TEXT PRIMARY KEY,
            experiment_id TEXT NOT NULL,
            pbs INTEGER NOT NULL,
            rtt INTEGER NOT NULL,
            FOREIGN KEY (experiment_id) REFERENCES experiments(id)
        )
        ''')
        
        self.db_connection.commit()

    def save_to_database(self):
        if not self.db_connection:
            messagebox.showwarning("Warning", "Please connect to a database first")
            return
            
        if not self.selected_cells:
            messagebox.showwarning("Warning", "No cells selected")
            return
            
        # Generate unique ID for experiment
        experiment_id = str(uuid.uuid4())
        
        try:
            cursor = self.db_connection.cursor()
            
            # Insert experiment
            cursor.execute('''
            INSERT INTO experiments (
                id, name, date, variant, chromosome, genomic_location, edit, gene_orientation
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                experiment_id,
                self.experiment_name.get(),
                datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                self.inputs["Variant"].get(),
                self.inputs["Chromosome"].get(),
                self.inputs["Genomic Location"].get(),
                self.inputs["Edit"].get(),
                self.inputs["Gene Orientation"].get()
            ))
            
            # Insert entries
            min_val = int(self.min_var.get())
            for row, col in self.selected_cells:
                entry_id = str(uuid.uuid4())
                cursor.execute('''
                INSERT INTO experiment_entries (
                    id, experiment_id, pbs, rtt
                ) VALUES (?, ?, ?, ?)
                ''', (
                    entry_id,
                    experiment_id,
                    min_val + col,  # PBS length
                    min_val + row   # RTT length
                ))
            
            self.db_connection.commit()
            messagebox.showinfo("Success", f"Experiment '{self.experiment_name.get()}' saved to database")
            
            # Generate a new unique experiment name for next save
            self.experiment_name.set(f"Experiment_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            
        except Exception as e:
            self.db_connection.rollback()
            messagebox.showerror("Error", f"Failed to save to database: {str(e)}")

def main():
    root = ttk.Window(themename="darkly", minsize=(1200, 645))  # Wider window, shorter height
    app = GridSelectorApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
