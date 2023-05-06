"""
Defines functions regarding outputting dataframes to latex
"""

#%%
#import csv
import subprocess
import os
import re
import math
import shutil
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

pd.set_option("display.max_colwidth", 1500)

#Color settings for figures (note: k = Black)
plt.rcParams.update({
        "figure.facecolor":     "white",
        "axes.facecolor":       "white",
        "axes.edgecolor":       "k",
        "savefig.facecolor":    "white",
        "axes.labelcolor":      "k",
        'legend.facecolor':     "white",
        "patch.edgecolor":      "k",
        "xtick.color":          "k",
        "ytick.color":          "k",
        "text.color":           "k",
})

#%% Functions
class TexTable:
    """The TexTable class is a dataframe coupled with captions, axis_labels, subplot titles, notes etc. 
    A list of TexTables can be used to print TexTables to a pdf using latex, either as tables or as figures."""
    def __init__(self, df, caption=None, note=None, label=None, decimals=2, use_index=False):
        self.data = df
        
        if caption is None:
            self.caption = "".join(e for e in df.columns[0] if e.isalnum()) + " against " + str(len(df.columns) - 1) + \
                           " y-variables: " + ", ".join(["".join(e for e in string if e.isalnum()) for string in list(df.columns[1:])])
        elif isinstance(caption, str):
            self.caption = caption
        #elif re.compile("[_]").search(caption) is None:
        #    self.caption = caption
        else:
            raise Exception("User-provided caption contains non-allowed special characters.")

        self.note = note
        self.decimals = decimals
        self.use_index = use_index


        if label is None:
            self.label = "temp_tbl"
        else:
            self.label = label.replace(" ", "_")


class TexFigure:
    """Class object for Figures intended to be LateX-output"""

    def __init__(self, fig, caption=None, note=None, label=None):

        if isinstance(fig, matplotlib.figure.Figure):
            self.fig = fig
        else:
            raise TypeError("The figure supplied is not of type matplotlib.figure.Figure.")
        
        self.note = note
        self.caption = caption
        
        if label is None:
            self.label = "temp_fig"
        elif isinstance(label, str):
            self.label = label.replace(" ", "_")

def initialize(filename):
    """A few checks to make sure the following code can run"""
    if not os.getcwd().split("\\")[-1] == "code":
        raise Exception("Current directory is not a 'code' folder")

    if any([False if letter.isalnum() or letter in [" ", "-", "_"] else True for letter in filename]):
        raise Exception("The chosen filename contains non-allowed special characters.")
        #replace("-", "a")
    # try:
    #     with open("../output/" + filename + ".pdf", "w+"):
    #         print("Printing to " + filename + ".pdf")
    # except PermissionError:
    #     raise PermissionError("Could not open file " + filename + ".pdf. Please close it before writing to it.")

    #Change settings such that thousand separator is ","
    pd.options.display.float_format = "{:,}".format

def compile_latex(filename):
    """ Compiles the .tex files to construct the pdfs."""

    #Run twice to make sure references are also updated in Table of contents
    first_run = subprocess.Popen(['pdflatex', filename + ".tex"], shell=False)
    first_run.wait()
    second_run = subprocess.Popen(['pdflatex', filename + ".tex"], shell=False)
    second_run.wait()

def clean_up(filename=False):
    """Removes unnecessary files from drive + moves pdfs"""
    if filename is not False:
        #Move PDF to output folder
        shutil.move(os.getcwd() + "\\" + filename + ".pdf", os.getcwd() + "\\..\\output\\" + filename + ".pdf")
    
    file_extensions = [".aux", ".log", ".lot", ".pdf", ".lof", ".tex", ".out"]

    #Delete .aux, .log etc. files after PDF has been generated:
    latex_files = [file for file in os.listdir(os.getcwd()) if file.endswith(tuple(file_extensions))]
    for file in latex_files:
        try:
            os.remove(os.getcwd() + "\\" + file)
        except OSError:
            print("The file" + file + " is still active and was not deleted.")


def insert_commas(num):

    """ Function used in write_figs_to_latex. It adds thousands separators to integers, because this cannot be done with regular formatting, 
    without keeping them as floats, in which case they have trailing ".0" which is ugly.""" 

    assert isinstance(num, str), "Number is not a string" 

    if ((len(num) > 2) and (num[0] == "-")):
        is_negative = True
        num = num[1:]
    else:
        is_negative = False

    split = num.split(".")
    assert len(split) == 1 or len(split) == 2, "Unexpected more than one decimal (period) in number"

    #The number before the decimal point
    num = split[0]

    try:
        int(num)
        is_number = True
    except ValueError:
        is_number = False
        
    if is_number:
        if len(num) > 4:
            n_separators = math.ceil(len(num)/3) - 1
            digit_groups = [num[i:i+3] for i in range(len(num)-n_separators*3, len(num), 3)]
            new_num = num[0:len(num)-n_separators*3] + "," + ",".join(digit_groups)
        else:
            new_num = num
        if len(split) == 2:
            new_num = new_num + "." + split[1]
        if is_negative:
            new_num = "-" + new_num
    elif num == "nan":
        new_num = "-"
    else:
        raise Exception("Number " + str(num) + " is not compatible")
    return new_num

def escape_underscores(string):
    """Escape underscores in regular text, keep if inside math mode where
    tex can actually compile them"""
    #Numbers should be returned as numbers
    if isinstance(string, str):
        #Inside math should not get subbed. The loop makes sure that multiple underscores in the same math mode get replaced.
        for _ in range(10):
            string = re.sub("(\\$.*?)_(.*?\\$)", "\\1!UNDERSCORE!\\2", string)
        #Inside \cref should net get subbed:
        string = re.sub("(\\\\cref{.*)_(.*})", "\\1!UNDERSCORE!\\2", string)
        #All other underscores are escaped
        string = re.sub("_", "\\_", string)
        #Bring back the underscores inside math mode
        string = string.replace("!UNDERSCORE!", "_") 
    return string

def latex_preamble(figures=False, tables=False):
    """ Set up the latex preamble (packages etc.)"""

    assert sum([figures, tables]) == 1, "Choose Figures or Tables"

    latex_start = "\n" + \
    "\\documentclass[12pt]{report}\n" + \
    "\\usepackage[colorlinks=true, breaklinks, linkcolor=black]{hyperref}\n" + \
    "\\usepackage[utf8]{inputenc}\n" + \
    "\\usepackage[T1]{fontenc}\n" + \
    "\\usepackage{threeparttable}\n" + \
    "\\usepackage{booktabs}\n" + \
    "\\usepackage{geometry}\n" + \
    "\\usepackage{floatrow}\n" + \
    "\\usepackage{amsmath}\n" + \
    "\\usepackage{amstext}\n" + \
    "\\usepackage{bm}\n" + \
    "\\floatsetup[table]{capposition=top}\n" + \
    "\\floatsetup[figure]{capposition=top}\n" + \
    "\\usepackage{float}\n" + \
    "\\usepackage{caption}\n" + \
    "\\captionsetup[table]{width=0.85\\textwidth}" + \
    "\\captionsetup[figure]{width=0.85\\textwidth}" + \
    "\\usepackage{graphicx}\n" + \
    "\\geometry{a4paper, top = 30mm, bottom = 30mm, left = 5mm, right = 5mm}\n" + \
    "\n" + \
    "\\usepackage{tocloft}\n" + \
    "\\renewcommand\\cftfigafterpnum{\\vskip7pt\\par}\n" + \
    "\\renewcommand\\cfttabafterpnum{\\vskip7pt\\par}\n" + \
    "\n" + \
    "\\usepackage[noabbrev, nameinlink]{cleveref}" + \
    "\n" + \
    "\\renewcommand*{\\arraystretch}{1.6}\n" + \
    "\n" + \
    "\\overfullrule=2cm\n" + \
    "\n" + \
    "\\begin{document}\n\n"

    latex_start = latex_start + \
    "All tables and figures comply with the data disclosure requirements of Statistics Denmark.\n" + \
    "No point or cell includes less than five observations, and no two observations are dominating.\n" + \
    "\\newpage\n"

    if figures:
        latex_start = latex_start + \
        "\\listoffigures\n" + \
        "\\thispagestyle{empty}\n" + \
        "\\newpage\n\n\n"
        
    elif tables:
        latex_start = latex_start + \
        "\\listoftables\n" + \
        "\\thispagestyle{empty}\n" + \
        "\\newpage\n\n\n"

    return latex_start

def write_figs_to_latex(list_of_TexFigures, filename):
    """
    Write dataframes to a single pdf as figures. Makes one figure for each TexTable, 
    and 1 subplot for each column in the TexTable (excluding the first, which is used as the x-axis)
    """

    #Add heading to PDF name
    filename = "FIGURES - " + filename.replace(".pdf", "")

    initialize(filename)
    
    latex_start = latex_preamble(figures=True, tables=False)

    latex_fig_start = "" + \
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" + \
    "%%%%%%%%%%%%%%%% START OF FIGURE %(label)s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" + \
    "\\begin{figure}[!htbp]\n" + \
    "\\centering\n" + \
    "\\includegraphics[scale=0.8]{%(label)s}\n" + \
    "\\caption{%(caption)s}\n" + \
    "\\label{fig:%(label)s}\n"

    latex_fig_end = "\\end{figure}\n\n" + \
    "%%%%%%%%%%%%%%%% END OF FIGURE %%%%%%%%%%%%%%%%%%%%\n\n\n\n\n\n"

    latex_end = "\\end{document}\n"

    
    with open(filename + ".tex", 'w') as f:
        f.write(latex_start)
        for fig_number, TexFig in enumerate(list_of_TexFigures):

            #Add figure number if it is a temporary figure
            if TexFig.label == "temp_fig":
                send_home = False
                prefix = str(fig_number) + "_"
            else:
                send_home = True
                prefix = ""
            
            #Save figure to drive
            TexFig.fig.savefig(prefix + TexFig.label + ".pdf", bbox_inches="tight", transparent=True, facecolor="none", edgecolor="none") #get.facecolor()
            
            if send_home:
                TexFig.fig.savefig("../output/send_home/" + prefix + TexFig.label + ".pdf", bbox_inches="tight", transparent=True)

            #Add note if it is written
            if TexFig.note is not None:
                latex_fig_start_fig = latex_fig_start + \
                "\\floatfoot{\\scriptsize \\emph{Notes:} " + escape_underscores(TexFig.note)  + "}\n"
                
            else:
                latex_fig_start_fig = latex_fig_start
            
            #Write everything to .tex file
            f.write(latex_fig_start_fig %{"caption": TexFig.caption, "fig_number": str(fig_number), "label": prefix + TexFig.label})
            f.write(latex_fig_end)

            if send_home:
                with open("../output/send_home/" + TexFig.label + ".tex", "w") as send_home_file:
                    send_home_file.write(latex_fig_start_fig %{"caption": TexFig.caption, "fig_number": str(fig_number), "label": prefix + TexFig.label})
                    send_home_file.write(latex_fig_end)

        f.write(latex_end)

    compile_latex(filename=filename)
    clean_up(filename)
    

def check_multicolumn_structure(multiindex):
    if (multiindex[0][0] != "") & (multiindex[0][0] == multiindex[1][0]):
        return "starts_with_a_multicolumn"
    else:
        return "starts_with_a_singlecolumn"

def write_tbls_to_latex(list_of_TexTables, filename):

    """ 
    Write dataframes to a single PDF. They must be TexTable objects and be gathered in a list.
    """
    #Add heading to PDF name
    filename = "TABLES - " + filename.replace(".pdf", "")

    initialize(filename)

    latex_start = latex_preamble(figures=False, tables=True)

    latex_end = "\\end{document}"

    latex_table_start = \
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" + \
    "%%%%%%%%%%%%%%%% START OF TABLE %(label)s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" + \
    "\\begin{table}[!htbp]\n" + \
    "\\centering\n" + \
    "\\caption{%(caption)s}\n" + \
    "\\label{tab:%(label)s}\n" + \
    "\\begin{threeparttable}\n" + \
    "\\scriptsize\n"

    

    latex_table_end_2 = "" + \
    "\\end{threeparttable}\n" + \
    "\\end{table}\n" + \
    "%%%%%%%%%%%%%%%% END OF TABLE %%%%%%%%%%%%%%%%%%%%\n\n\n\n\n\n"
    

    with open(filename + ".tex", 'w') as f:
        # nan_val = "-"
        f.write(latex_start)
        for tbl_number, tbl in enumerate(list_of_TexTables):
            if tbl.label == "temp_tbl":
                send_home = False
                prefix = str(tbl_number) + "_"
            else:
                send_home = True
                prefix = ""
            
            #Write table setup
            f.write(latex_table_start %{"caption": tbl.caption, "label": prefix + tbl.label})

            #Format content of table
            df = tbl.data.copy()
            for col in df.columns:

                #Identify rows with ints, floats and strings respectively 
                int_rows = df.loc[:, col][(df.loc[:, col].map(type) == int)].index
                text_rows = df.loc[:, col][(df.loc[:, col].map(type) == str)].index
                float_rows = df.loc[:, col][(df.loc[:, col].apply(lambda x: isinstance(x, float)))].index

                #Insert commas in ints above 4 digits
                df.loc[int_rows, col] = df.loc[int_rows, col].astype("str").apply(insert_commas)

                #Format floats
                if df.loc[float_rows, col][df.loc[float_rows, col].notnull()].apply(float.is_integer).all():
                    df.loc[float_rows, col] = df.loc[float_rows, col].astype("Int64").astype("str").apply(insert_commas)
                else:
                    df.loc[float_rows, col] = df.loc[float_rows, col].astype(float).round(tbl.decimals).apply(("{0:." + str(tbl.decimals) + "f}").format).apply(insert_commas)
                
                #apply(lambda x: str("{:,.2f}".format(round(x, 2))) if not math.isnan(x) else nan_val)

                #Text rows: escape underscores (escapes outside of math, keeps inside)
                df.loc[text_rows, col] = df.loc[text_rows, col].apply(escape_underscores)

            #Escape special characters in column names where appropriate
            if not isinstance(df.columns, pd.MultiIndex):
                for _, name in enumerate(df.columns):
                    if isinstance(name, str):
                        newname = escape_underscores(name)
                        newname = re.sub("(%|&)", "backslash\\1", newname).replace("backslash", "\\")
                        df.rename({name:newname}, axis="columns", inplace=True)
            else:
                cmidlines = get_cmidlines(df.columns)
                for level in range(0, df.columns.nlevels):
                    #A bug in pd.DataFrame.to_latex means it cannot handle empty column labels
                    d = dict(zip(df.columns.levels[level], ["!PLACEHOLDER!" if c in ["", " ", "   ", "_"] else c for c in df.columns.levels[level]]))
                    df = df.rename(columns=d, level=level)
                    
            #Remove index names, since they can ruin the LateX output (midrules in wrong spots)
            if isinstance(df.index, pd.Index):
                df.index.name = None
            else:
                df.index.names = [None]*len(df.index.levels)
            
            #Write table content
            latex_table = df.to_latex(index=tbl.use_index, column_format="l"*(len(df.columns)+tbl.use_index), multicolumn_format="c", escape=False)
            # Ændre her så index = true kan bruges.
                # + index ves column_format.
                # Brug det og skriv et ekstra l i begin{tabular}

            if isinstance(df.columns, pd.MultiIndex):
                if check_multicolumn_structure(df.columns) == "starts_with_a_multicolumn":
                    line = latex_table.split("\n")[2]
                    pattern = r"multicolumn{([0-9])}"
                    missing_cols = len(df.columns) - sum([int(t) for t in re.findall(pattern, line)])
                    insert = r"\\multicolumn{" + str(missing_cols) + r"}{c}{\1} \2"
                    if "&" not in line:
                        pattern = r"(^.*?) (\\\\)"
                        # insert = insert + r"\\\\"
                    elif "&" in line:
                        pattern = r"(^.*?) (&)"
                        # insert = insert + r"&"
                    latex_table = latex_table.replace(line, re.sub(pattern, insert, line, count=1))
                # elif check_multicolumn_structure(df.columns) == "starts_with_a_singlecolumn":
                latex_table = latex_table.replace("!PLACEHOLDER!", "")

            #Make multicolumn where first column says "Panel A/B/C/..." across all columns if they are empty
            insert = r"\\" + "multicolumn{" + str(len(df.columns)) + r"}{l}{\\textbf{\1}}            "
            pattern = "(Panel [A-Z]{1}[^&]+)(&)([ ]+)" + "(&[ ]+)"*(len(df.columns) - 2)
            latex_table = re.sub(pattern, insert, latex_table)

            #Insert midrules in rows with MIDRULE in first column (and rest empty)
            insert = r"\\" + "midrule"
            pattern = "MIDRULE[^&]+&[ ]+" + "&[ ]+"*(len(df.columns) - 2) + r"\\\\"
            latex_table = re.sub(pattern, insert, latex_table)

            # Insert midrule over "Observations" row
            # Insert dobblet rule over første række der starter med "F-stat" og lige over søjlerne? 

            if isinstance(df.columns, pd.MultiIndex):
                #Make lines underneath categories/multicolumns
                line_pattern = r"toprule\n.*\\\\\n"
                match = re.search(line_pattern, latex_table).group(0)
                latex_table = latex_table.replace(match, match + cmidlines)

            f.write(latex_table)

            if tbl.note is not None:
                latex_table_end_tbl = "\\begin{tablenotes}[flushleft]\n\\scriptsize\n" + "\\item \\emph{Notes:}\n" + \
                                    escape_underscores(tbl.note) + " \n\\end{tablenotes}\n" + latex_table_end_2
            else:
                latex_table_end_tbl = latex_table_end_2
            f.write(latex_table_end_tbl)

            #Write the same things if we send home (if label is defined)
            if send_home:
                with open(os.getcwd() + "/../output/send_home/" + tbl.label + ".tex", "w") as send_home_file:
                    send_home_file.write(latex_table_start %{"caption": tbl.caption, "label": prefix + tbl.label})
                    send_home_file.write(latex_table)
                    send_home_file.write(latex_table_end_tbl)

        f.write(latex_end)

    compile_latex(filename=filename)
    clean_up(filename)

def get_cmidlines(columns):
    """ Construct the line in latex of the form "\cmidrule(lr){2-4}\cmidrule..." which implements lines under grouped columns"""
    assert isinstance(columns, pd.MultiIndex)
    res = ""
    at_start = True
    curr_title = None
    col_count = 1
    repetition_count = 1
    for i, pair in enumerate(columns.to_list()):
        # print(i, pair[0])
        # print(col_count)
        if (at_start) & ((pair[0] == "") or (pair[1] == "")):
            # res += "empty_"
            col_count += 1
        else:
            if pair[0] == "_":
                col_count += 1
                continue
            if curr_title == pair[0]:
                repetition_count += 1            
            else:
                curr_title = pair[0]
                if at_start:
                    at_start = False
                else:
                    res += "\\cmidrule(lr){" + str(col_count) + "-" + str(col_count + repetition_count - 1) + "}"
                    col_count += repetition_count
                repetition_count = 1
            #if last index, end it
            if i == len(columns) - 1:
                res += "\\cmidrule(lr){" + str(col_count) + "-" + str(col_count + repetition_count - 1) + "}"
    return res + "\n"


#%%
