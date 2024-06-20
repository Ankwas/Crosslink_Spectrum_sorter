# ----- Packages ----- 
import customtkinter as ctk
from tkinter import filedialog
from tkinter import messagebox
import numpy as np
from pyteomics import mgf
# import subprocess

# ----- Functions ----------
def find_pairs(tolerance ,mz_values, intensity_values, target_difference, intensity_cutoff, low_mass_cutoff):
    """
    Finds pairs of m/z values with a specified difference and intensity cutoff.
    Args:
        tolerance (float): The tolerance on the intenseties
        mz_values (numpy.array): Array of m/z values.
        intensity_values (numpy.array): Array of intensity values.
        target_difference (float): Desired difference between m/z values.
        intensity_cutoff (float): Intensity cutoff threshold.
        low_mass_cutoff (float): Low mass cutoff threshold.  
    Returns:
        list: List of tuples representing pairs of m/z values and intensities.
    """
    pairs = []
    # Apply intensity cutoff
    valid_intensity = (intensity_values >= intensity_cutoff * intensity_values.max())
    mz_values_filtered = mz_values[valid_intensity]
    intensity_values_filtered = intensity_values[valid_intensity]
    # Apply low mass cutoff
    low_mass_indices = mz_values_filtered >= low_mass_cutoff
    mz_values_filtered = mz_values_filtered[low_mass_indices]
    intensity_values_filtered = intensity_values_filtered[low_mass_indices]
    #print("Filtered mz_values:", mz_values_filtered)  # Debug print
    sorted_indices = np.argsort(mz_values_filtered)
    mz_values_filtered = mz_values_filtered[sorted_indices]
    intensity_values_filtered = intensity_values_filtered[sorted_indices]
    for i, mz1 in enumerate(mz_values_filtered):
        for j, mz2 in enumerate(mz_values_filtered[i+1:]):
            diff = mz2 - mz1
            # print("Difference:", diff)  # Debug print
            if abs(diff - target_difference) <= tolerance*2:
                intensity1 = intensity_values_filtered[i]
                intensity2 = intensity_values_filtered[i+j+1]
                pairs.append((mz1, intensity1, mz2, intensity2))
    # print("Found pairs:", pairs)  # Debug print
    return pairs

def write_spectra_with_pairs(writer, spectrum, pairs):
    """
    Writes spectra with pairs to the output MGF file.
    Args:
        writer (file): File object for writing.
        spectrum (dict): Spectrum data.
        pairs (list): List of pairs of m/z values and intensities.
    """
    writer.write("BEGIN IONS\n")
    for key, value in spectrum['params'].items():
        if isinstance(value, tuple):
            writer.write('{}={} {}\n'.format(key.upper(), value[0], value[1]))
        else:
            writer.write('{}={}\n'.format(key.upper(), value))
    for mz1, intensity1, mz2, intensity2 in pairs:
        writer.write( "{:.4f} {:.4f}\n".format(mz1, intensity1)) # {:.4f} is to print with only 4 significant digits
        writer.write( "{:.4f} {:.4f}\n".format(mz2, intensity2))
    writer.write("END IONS\n")

def main(entry_target_difference, entry_intensity_cutoff, entry_low_mass_cutoff, entry_min_charge,
         entry_min_pair, entry_tolerance, entry_input_mgf_file, entry_output_mgf_file):
    """
    Runs the script to process input MGF files and extract spectra with desired properties.
    Args:
        entry_target_difference (float), value for the distance between to peaks in the spectrum for them to be a pair
        entry_intensity_cutoff (float), cutoff value for intensetys, number is in procent 
        entry_low_mass_cutoff (float), cutoff value for mass, number is in m/z 
        entry_min_charge (float), minimum charge a peptide needs to be taken into acount
        entry_min_pair (float), The minimum number of pairs in the peptide 
        entry_tolerance (float), the tolerance of the intensety values
        entry_input_mgf_file (string), locations from where the input file is loaded
        entry_output_mgf_file (string), location to where the output file is saved
    """
    # Load data from the UI 
    try:
        target_difference = float(entry_target_difference)
    except ValueError:
        messagebox.showerror("Error", "Please fill in taget difference with valid numbers.")
    try: 
        intensity_cutoff = float(entry_intensity_cutoff)/100  # Divide by 100 to get the intensity as 0.1 instead of 10
    except ValueError:
        messagebox.showerror("Error", "Please fill in taget intensety cutoff with valid numbers.")
    try:
        low_mass_cutoff = float(entry_low_mass_cutoff)
    except ValueError:
        messagebox.showerror("Error", "Please fill in low mass cutoff with valid numbers.")
    try:
        min_charge = float(entry_min_charge)
    except ValueError:
        messagebox.showerror("Error", "Please fill in minimum charge with valid numbers.")
    try:
        min_pair = float(entry_min_pair)
    except ValueError:
        messagebox.showerror("Error", "Please fill in minumum nuber of pairs with valid numbers.")
    try:
        tolerance = float(entry_tolerance)
    except ValueError:
        messagebox.showerror("Error", "Please fill in tolerance.")
    input_mgf_file = entry_input_mgf_file
    output_mgf_file = entry_output_mgf_file
    # Check if data is loaded 
    if not (input_mgf_file and output_mgf_file):
        messagebox.showerror("Error", "Please select input and output files.")
        return
    # Process input MGF files 
    try:
        with mgf.read(input_mgf_file) as reader:
            with open(output_mgf_file, 'w') as writer:
                for spectrum in reader:
                    charge = int(spectrum['params']['charge'][0])
                    if charge < min_charge:
                            continue
                    mz_values = np.array(spectrum.get('m/z array', []))
                    intensity_values = np.array(spectrum.get('intensity array', []))
                    mz_int_dict = dict(zip(mz_values, intensity_values))
                    pairs = find_pairs(tolerance, mz_values, intensity_values, target_difference, intensity_cutoff, low_mass_cutoff)
                    if len(pairs) > min_pair:
                        print(pairs, "\n", "=" * 10) #debug print 
                        write_spectra_with_pairs(writer, spectrum, pairs)
        messagebox.showinfo("Success", "Script executed successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")

def ui():
    """
    Creates the userinterface that the user will se when they run the program
    Contains sub functions that alows for user inputs and initiates the program
    """

    # Subfunctions for the userinterface
    def browse_input_mgf_file():
        """
        Makes a button that opens a file explore that only shows mgf files to select as input
        """
        filename = filedialog.askopenfilename(filetypes=[("MGF files", "*.mgf")])
        entry_input_mgf_file.delete(0, ctk.END)
        entry_input_mgf_file.insert(0, filename)

    def browse_output_mgf_file():
        """
        Makes a button that opens a file explorer to select the output file destination.
        The output file name will be the input file name + '_sorted.mgf'.
        """
        input_filename = entry_input_mgf_file.get()
        output_filename = filedialog.asksaveasfilename(initialfile=input_filename.split('/')[-1].replace('.mgf', '_sorted.mgf'), defaultextension=".mgf", filetypes=[("MGF files", "*.mgf")])
        entry_output_mgf_file.delete(0, ctk.END)
        entry_output_mgf_file.insert(0, output_filename)

    def set_target_difference():
        """
        Makes a button that sets the input to the distance from one BuUrBu peak to another
        """
        entry_target_difference.delete(0, ctk.END)
        entry_target_difference.insert(0, "25.97")

    def set_min_charge():
        """
        Makes a button that sets the input to 3
        """
        entry_min_charge.delete(0, ctk.END)
        entry_min_charge.insert(0, "3")

    def run_program():
        """
        takes all the inputs from the user interface and gives to the program
        """
        main(entry_target_difference.get(),
             entry_intensity_cutoff.get(),
             entry_low_mass_cutoff.get(),
             entry_min_charge.get(),
             entry_min_pair.get(),
             entry_tolerance.get(),
             entry_input_mgf_file.get(),
             entry_output_mgf_file.get())
    
    # User interface setup
    ctk.set_appearance_mode("light")
    ctk.set_default_color_theme("dark-blue")

    # Create main window
    root = ctk.CTk()
    root.title("Spektrum sorter")

    # Create input fields and buttons
    label_target_difference = ctk.CTkLabel(root, text="Target Difference:")
    label_target_difference.grid(row=0, column=0)
    entry_target_difference = ctk.CTkEntry(root)
    entry_target_difference.grid(row=0, column=1)
    button_set_target_difference = ctk.CTkButton(root, text="Set to BuUrBu", command=set_target_difference)
    button_set_target_difference.grid(row=0, column=2)

    label_intensity_cutoff = ctk.CTkLabel(root, text="  Intensity Cutoff (ex. 10 for cut of at 10% intensety):")
    label_intensity_cutoff.grid(row=1, column=0)
    entry_intensity_cutoff = ctk.CTkEntry(root)
    entry_intensity_cutoff.grid(row=1, column=1)

    label_low_mass_cutoff = ctk.CTkLabel(root, text="Low mass cutoff (ex. 400):")
    label_low_mass_cutoff.grid(row=2, column=0)
    entry_low_mass_cutoff = ctk.CTkEntry(root)
    entry_low_mass_cutoff.grid(row=2, column=1)

    label_min_charge = ctk.CTkLabel(root, text="Minimum charge:")
    label_min_charge.grid(row=3, column=0)
    entry_min_charge = ctk.CTkEntry(root)
    entry_min_charge.grid(row=3, column=1)
    # Button for the trail and error part of the program where layzenes made this easier than type 3 
    # button_set_min_charge = ctk.CTkButton(root, text="Set to 3", command=set_min_charge)
    # button_set_min_charge.grid(row=3, column=2)

    label_min_pair = ctk.CTkLabel(root, text="Minimum number of pairs:")
    label_min_pair.grid(row=4, column=0)
    entry_min_pair = ctk.CTkEntry(root)
    entry_min_pair.grid(row=4, column=1)

    label_tolerance = ctk.CTkLabel(root, text="Tolerance of spectrum peaks in Da (ex 0.1)")
    label_tolerance.grid(row=5, column=0)
    entry_tolerance = ctk.CTkEntry(root)
    entry_tolerance.grid(row=5, column=1)

    label_input_mgf_file = ctk.CTkLabel(root, text="Input MGF File:")
    label_input_mgf_file.grid(row=6, column=0)
    entry_input_mgf_file = ctk.CTkEntry(root)
    entry_input_mgf_file.grid(row=6, column=1)
    button_browse_input_mgf = ctk.CTkButton(root, text="Browse", command=browse_input_mgf_file)
    button_browse_input_mgf.grid(row=6, column=2)

    label_output_mgf_file = ctk.CTkLabel(root, text="Output MGF File:")
    label_output_mgf_file.grid(row=7, column=0)
    entry_output_mgf_file = ctk.CTkEntry(root)
    entry_output_mgf_file.grid(row=7, column=1)
    button_browse_output_mgf = ctk.CTkButton(root, text="Browse", command=browse_output_mgf_file)
    button_browse_output_mgf.grid(row=7, column=2)

    button_run = ctk.CTkButton(root, text="Run Script",command=run_program)
    button_run.grid(row=8, column=1)

    # Starts the userinterface
    root.mainloop()

ui()
