# SCP_LF
Single cell workflow for label free with Jupyter notebook and data manage system

# How to Use
## Installing Github Jupiter Notebook Workflow:  
1.	Download Visual Studio Code with Python  
2.	Download and unzip the repository from https://github.com/RTKlab-BYU/SCP_LF  
3.	Open the folder in VS Code  
4.	Open a new terminal and type “pip install -r requirements.txt’ and then hit enter.  
## Analyzing the Proteome:  
1.	Place the output files from Proteome Discoverer in your unzipped folder. make sure Abundances was exported from PD    
2.	Open Single_cell_workflow.ipynb.  
3.	First decide if you want automatically exported files, if so near the top of the notebook change WRITE_OUTPUT = True.  
a.	Average ID’s and the corresponding standard deviations and 95% confidence intervals, CVs, protein lists, principal components, and volcano graph comparisons will be exported.  
4.	Name your input files:  
a.	Change the dictionary in the filelist variable which is near the bottom of the notebook to have the names of your Proteome Discoverer output.  
i.	“input1” is your Proteins.txt file or combined_proteins.tsv or diann-output.pg_matrix.tsv  
ii.	“input2” is your PeptideGroups.txt file or combined_peptides.tsv or diann-output.pr_matrix.tsv  
i.	“input3” is anything for PD or combined_proteins.tsv or diann-output.pg_matrix.tsv  
ii.	“input4” is anything for PD or combined_peptides.tsv without MBR or diann-output.pr_matrix.tsv without MBR.  
iv.	“input5” is your InputFiles.txt file.  
v.	“process_app” should be set to “PD” for Proteome Discoverer, "FragPipe" for Fragpipe or "DIANN" for DIANN.  
b.	To add multiple analyses add more dictionaries in curly brackets separated by commas with the same variables defined.  
c.	Change SETTINGS_FILE to the name of the file   
d.	If you place any PD files in subdirectories please add the sub directory followed by “/”  
5.	Set the SETTINGS_FILE variable to be equal to the file name you make in the next step enclosed by quotation marks.  
6.	Create the SETTINGS_FILE  
   a.	In Excel or another spreadsheet application open “settings_test.txt”.  
   b.	The “Conditions” should contain unique identifiers that will be used in most graphs.  
   c.	The filter_in column will be equal to the filename convention you used in naming your raw files    
   i.	If you used multiple Proteome Discoverer analyses in step 4b, then add @ and the index of the analysis containing the files for this group starting with 0, not 1.    
   d.	filter_out should be equal to the pattern of any files you don’t want included in this group that match the filter_in convention, if any. Can be left blank if there are none.    
   e.	Any other columns are optional, unless you are using grouped, stacked, or grouped_stacked ID/CV graphs, with any grouping besides “ID_Mode”. In that case column names can be set as the grouping variables they will appear in the order they are in the SETTINGS FILE  
   f.	Save as a tab separated txt within your unzipped folder    
8.	Configure Graphical Output at the bottom of the notebook and rerun these blocks individually until you have the desired output.  
   a.	ID Plots  
        i.	Specify if mean labels should be shown (True or False).  
        ii.	Specify whether to hide error bars (None) or use confidence interval (ci95), or standard deviation for error bars (stdev).  
        iii.	Change X and Y titles to match your data (make sure to change Protein to peptide for peptide graphs on the Y axis.  
        iv.	Specify colors to be used using either hex codes or color names.  
        v.	Specify graph width and height.  
        vi.	Specify font.  
        vii.	Specify graph type.  
        viii.	Specify grouping variables (these are the column names or “ID_Mode” to use MS2 vs MBR identifications.  
        ix.	Specify plot type (1 for protein, 2 for peptide).  
   b.	CV Plots  
        i.	Specify if mean labels should be shown (True or False).  
        ii.	Specify if there should be a boxplot within the violin plot  
        iii.	Change X and Y titles to match your data (make sure to change Protein to peptide for peptide graphs on the Y axis.  
        iv.	Specify colors to be used using either hex codes or color names.  
        v.	Specify graph width and height.  
        vi.	Specify font.  
        vii.	Specify y axis range.  
        viii.	Specify graph type.  
        ix.	Specify grouping variables (these are the column names or “ID_Mode” to use MS2 vs MBR identifications.  
        x.	Specify plot type (1 for protein, 2 for peptide).  
   c.	ID venn plot  
        i.	Specify groups to compare (up to 3), using the values in the Group Name column in the SETTINGS_FILE.  
        ii.	Input chart title.  
        iii.	Change opacity of circles if desired.  
        iv.	Specify colors for groups.  
   d.	Volcano Plot  
        i.	Specify groups to compare (up to 2), using the values in the Group Name column in the SETTINGS_FILE.  
        ii.	Input chart title.  
        iii.	Specify axis titles.  
        iv.	Specify color for significantly up-regulated proteins.  
        v.	Specify color for significantly down-regulated proteins.  
        vi.	Specify color for proteins without significant enrichment between samples.  
        vii.	Specify plot height and width.  
        viii.	Specify font.  
        ix.	Specify axis ranges if desired (perhaps you have outliers)  
   e.	Principal Component Analysis (PCA)  
        i.	 specify groups to compare (as many as you want), using the values in the Group Name column in the SETTINGS_FILE.  
        ii.	Input chart title.  
        iii.	Specify colors to be used using either hex codes or color names.  
        iv.	Specify symbols for each group.  
        v.	Specify marker size.  
        vi.	Specify graph height and width.  
        vii.	Specify font.  
        viii.	Specify axis ranges if desired (perhaps you have outliers)  
   f. Heatmap  
        i. log2 applies to view only  
    	  ii. significant calculated both ways for output file using ANOVA  
## Proteome Discoverer 2.5 Analysis and Export for Chromatography:
1.	Download apQuant from https://ms.imp.ac.at/?action=apQuant and install.
2.	Restart and open Proteome Discoverer 2.5 
3.	Open a new study
4.	Load the analysis by closing any workflows in your PD study and then clicking the “Analysis from Template” and load the analysis template provided in this paper.
5.	Ensure that you put the contaminants database in the protein marker node and the protein database into the protein annotation, Sequest HT, INFERYS, and spectrum files RC nodes. Notice there are two Sequest HT nodes.
6.	Drag the input files over to the analysis panel on the right.
7.	Start analysis and ignore warnings involving grouping.
8.	Make sure the “Identified By”, “Spectrum File”, ”Sequence”, “Peak Apex”, “FWHM”, “Apex Intensity” and “Area” columns are visible in the apQuant Features tab.
9.	Export the apQuant Features to a tab delimited txt file. (File>>Export>>To Text (tab delimited…)
10.	Files will be in the same folder as your study under a subfolder with your study’s name.
## Analyzing the Chromatography:
1.	Copy output of PD 2.5 with apQuant to your unzipped folder
2.	Open “Chromatography_workflow.ipynb”
3.	Near the bottom of the page, set PEAKS_DATA_FILES to be equal to be your filenames from PD 2.5 with apQuant, enclosed in quotation marks and separated by commas.
4.	Set the PEAKS_SETTINGS_FILE variable to be equal to the file name you make in the next step enclosed by quotation marks.
5.	Create the SETTINGS_FILE
a.	In Excel or another spreadsheet application open “settings_test.txt”.
b.	The “Group Name” should contain unique identifiers that will be used in most graphs.
c.	The filter_in column will be equal to the filename convention you used in naming your raw files.
i.	If you used multiple Proteome Discoverer analyses in step 4b, then add @ and the index of the analysis containing the files for this group starting with 0, not 1.
d.	filter_out should be equal to the pattern of any files you don’t want included in this group that match the filter_in convention, if any. Can be left blank if there are none.
e.	Any other columns are optional, unless you are using grouped graphs. In that case column names can be set as the grouping variables
6.	Configure plots:
a.	Decide if you want a median label in each box plot.
b.	Specify X and Y axis labels.
c.	Specify colors to be used.
d.	Specify if log10 transformations will be used.
e.	Specify width and height of graph.
f.	Specify Font.
g.	Specify whether outliers are removed:
i.	“remove” removes all outliers.
ii.	“suspectedoutliers” removes suspected outliers.
iii.	False leaves them in.
iv.	“all” shows all points.
v.	“outliers” removes all outliers.
h.	Specify axis y limits if you want to focus on part of the graph
i.	Specify mode and grouping to group by a column specified in SETTINGS_FILE.

