import wx
import wx.grid
import pandas as pd
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import io
import os
import glob
import csv

# Load data
annotation1 = pd.read_csv("source_data/annotation.txt",sep=",", low_memory=False)
annotation1 = annotation1.rename(columns={annotation1.columns[8]: "gene"})
full_list = pd.read_csv("source_data/full_list.txt",sep="\t", low_memory=False)
diseases = full_list['disease'].unique()
disease_type = full_list['disease_type'].unique()
gene = full_list['gene'].unique()

################ Tab 1: Biomarker ####################
######################################################
######################################################

class TabOne(wx.Panel):
        def __init__(self, parent):
            wx.Panel.__init__(self, parent)
            
            # Create panels
            side_panel = wx.Panel(self)
            side_panel.SetBackgroundColour('#D3D3D3')
            main_panel = wx.Panel(self)
            about_gene_panel = wx.Panel(main_panel)
            about_gene_panel.SetBackgroundColour('#FFFFFF')
            transcript_panel = wx.Panel(main_panel)
            transcript_panel.SetBackgroundColour("#FFFFFF")
            side_sizer = wx.BoxSizer(wx.VERTICAL)
            about_gene_sizer = wx.BoxSizer(wx.VERTICAL)
            transcript_sizer = wx.BoxSizer(wx.VERTICAL)
            
            # Disease selection
            disease_sizer = wx.BoxSizer(wx.VERTICAL)
            disease_label = wx.StaticText(side_panel, label='Select the disease of interest:')
            self.disease_combo = wx.ComboBox(side_panel, choices=list(diseases), style=wx.CB_READONLY)
            self.disease_combo.Bind(wx.EVT_COMBOBOX, self.OnDiseaseSelected)
            disease_sizer.Add(disease_label, 0, wx.ALL | wx.RIGHT, 5)
            disease_sizer.Add(self.disease_combo, 0, wx.ALL, 5)
            side_sizer.Add(disease_sizer, 0, wx.EXPAND)

            # Disease type selection
            self.disease_type_sizer = wx.BoxSizer(wx.VERTICAL)
            disease_type_label = wx.StaticText(side_panel, label='Select the disease type:')
            self.disease_type_combo = wx.ComboBox(side_panel, choices=list(disease_type), style=wx.CB_READONLY|wx.TE_MULTILINE)
            self.disease_type_combo.SetMinSize((250, -1))
            self.disease_type_combo.Bind(wx.EVT_COMBOBOX, self.OnDiseaseTypeSelected)
            self.disease_type_sizer.Add(disease_type_label, 0, wx.ALL | wx.RIGHT, 5)
            self.disease_type_sizer.Add(self.disease_type_combo, 0, wx.ALL, 5)
            side_sizer.Add(self.disease_type_sizer, 0, wx.EXPAND)

            # gene selection
            self.gene_sizer = wx.BoxSizer(wx.VERTICAL)
            gene_label = wx.StaticText(side_panel, label='Gene of Interest:')
            self.gene_combo = wx.ComboBox(side_panel, choices=list(gene), style=wx.CB_READONLY)
            self.gene_combo.Bind(wx.EVT_COMBOBOX, self.process_data)
            self.gene_sizer.Add(gene_label, 0, wx.ALL | wx.RIGHT, 5)
            self.gene_sizer.Add(self.gene_combo, 0, wx.ALL, 5)
            side_sizer.Add(self.gene_sizer, 0, wx.EXPAND | wx.ALL, 10)
            side_panel.SetSizerAndFit(side_sizer)
            
            # create table gene info 
            about_gene_label = wx.StaticText(about_gene_panel, label="About the gene")
            gene_label_font = about_gene_label.GetFont()
            gene_label_font.SetPointSize(16)
            about_gene_label.SetFont(gene_label_font)
            self.gene_info_table = wx.grid.Grid(about_gene_panel)
            self.gene_info_table.CreateGrid(numRows=1, numCols=6)
            self.gene_info_table.SetColLabelValue(0, "Gene")
            self.gene_info_table.SetColLabelValue(1, "MOI")
            self.gene_info_table.SetColLabelValue(2, "SOP")
            self.gene_info_table.SetColLabelValue(3, "Classification")
            self.gene_info_table.SetColLabelValue(4, "Online report")
            self.gene_info_table.SetColLabelValue(5, "Classification date")
            self.gene_info_table.SetRowSize(0, 50)
            about_gene_sizer.Add(about_gene_label, 0, wx.ALL | wx.RIGHT, 5)
            about_gene_sizer.Add(self.gene_info_table, 1, wx.ALL | wx.RIGHT, 5)
            self.gene_info_table.AutoSizeColumns()
            about_gene_panel.SetSizerAndFit(about_gene_sizer)
            
            # create table gene transcripts
            gene_trans_label = wx.StaticText(transcript_panel, label="Gene's transcript")
            trans_label_font = gene_trans_label.GetFont()
            trans_label_font.SetPointSize(16)
            gene_trans_label.SetFont(trans_label_font)
            self.gene_trans_table = wx.grid.Grid(transcript_panel)
            self.gene_trans_table.CreateGrid(numRows= 0, numCols=5)
            self.gene_trans_table.SetColLabelValue(0, "Gene")
            self.gene_trans_table.SetColLabelValue(1, "Transcript name")
            self.gene_trans_table.SetColLabelValue(2, "Transcript type")
            self.gene_trans_table.SetColLabelValue(3, "Transcription start site (TSS)")
            self.gene_trans_table.SetColLabelValue(4, "Transcript end (bp)")
            self.gene_trans_table.SetColLabelValue(5, "Transcript start (bp)")
            transcript_sizer.Add(gene_trans_label, 0, wx.ALL | wx.RIGHT, 5)
            transcript_sizer.Add(self.gene_trans_table, 4, wx.ALL | wx.RIGHT, 5)
            self.gene_trans_table.AutoSizeColumns()
            transcript_panel.SetSizerAndFit(transcript_sizer)
            
            
            # Main sizer to hold all panels
            main_sizer = wx.BoxSizer(wx.VERTICAL)
            main_sizer.Add(about_gene_panel, 0, wx.EXPAND)
            main_sizer.Add(transcript_panel, 1, wx.EXPAND)
            main_panel.SetSizerAndFit(main_sizer)
            
            sizer = wx.BoxSizer(wx.HORIZONTAL)
            sizer.Add(side_panel, 1, wx.EXPAND)
            sizer.Add(main_panel, 4, wx.EXPAND)
            self.SetSizer(sizer)
        
            

        def OnDiseaseSelected(self, event):
            # Update disease type
            selected_disease = self.disease_combo.GetValue()
            disease_types = full_list.loc[full_list['disease'] == selected_disease, 'disease_type'].unique()
            self.disease_type_combo.SetItems(list(disease_types))
            self.disease_type_sizer.ShowItems(show=True)
            
        def OnDiseaseTypeSelected(self, event):
            # Update gene list
            selected_disease_type = self.disease_type_combo.GetValue()
            gene_list = full_list.loc[full_list['disease_type'] == selected_disease_type, 'gene'].unique()
            self.gene_combo.SetItems(list(gene_list))
            self.gene_sizer.ShowItems(show=True)
        
            self.Layout()
            
        def process_data(self, event):
            # filter
            disease_n = self.disease_combo.GetValue()
            disease_t = self.disease_type_combo.GetValue()
            target_gene = self.gene_combo.GetValue()
            gene_info = full_list.loc[(full_list['disease'] == disease_n) & (full_list['disease_type'] == disease_t) &
                                    (full_list['gene'] == target_gene)]
            
            gene_info = gene_info.drop(columns=["GENE ID (HGNC)","disease_type","DISEASE ID (MONDO)","disease"])
            
            # Update the table with filtered data
            # Clear the existing grid data
            self.gene_info_table.ClearGrid()
            # Get the number of rows and columns from the grid
            num_rows, num_cols = gene_info.shape
            
            # Populate the grid
            for row in range(num_rows):
                for col in range(num_cols):
                    value = str(gene_info.iloc[row, col])
                    self.gene_info_table.SetCellValue(row, col, value)
            
            # Refresh the grid
            self.gene_info_table.Refresh()
            self.gene_info_table.AutoSizeColumns()
            
            # annotate genes
            gene_trans = annotation1[annotation1["gene"] == target_gene]
            gene_trans = gene_trans.iloc[:, [8, 9, 12, 13, 14, 15]]
            gene_trans = gene_trans.drop_duplicates()
            
            # Clear existing grid data
            self.gene_trans_table.ClearGrid()
            
            # Get the number of rows and columns from the data_frame
            num_rows, num_cols = gene_trans.shape
            # Get the current number of rows in the grid
            current_num_rows = self.gene_trans_table.GetNumberRows()
            current_num_cols = self.gene_trans_table.GetNumberCols()
            # Append additional rows to the grid
            if current_num_rows < num_rows:
                self.gene_trans_table.AppendRows(numRows=num_rows - current_num_rows)

            if current_num_cols < num_cols:
                self.gene_trans_table.AppendCols(numCols=num_cols - current_num_cols)

            # Populate the grid
            for row in range(num_rows):
                for col in range(num_cols):
                    value = str(gene_trans.iloc[row, col])
                    self.gene_trans_table.SetCellValue(row, col, value)
                    
            # Refresh the grid
            self.gene_trans_table.Refresh()
            self.gene_trans_table.AutoSizeColumns()
            
            self.Layout()
        
#####################################################################
######Tab 2: Expression##############################################
#####################################################################

class TabTwo(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        #self.current_page = 1
        # Create panels
        side_panel = wx.Panel(self)
        side_panel.SetBackgroundColour('#D3D3D3')
        main_panel = wx.Panel(self)
        side_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        #filter by biomarker by disease
        diseases = full_list['disease'].unique()
        disease_input_label = wx.StaticText(side_panel, label='Select the disease of interest:')
        self.disease_combo = wx.ComboBox(side_panel, choices=list(diseases), style=wx.CB_READONLY)

        # File input
        file_label = wx.StaticText(side_panel, label="Choose file to upload")
        self.file_input = wx.FilePickerCtrl(side_panel, style=wx.FLP_DEFAULT_STYLE | wx.FLP_USE_TEXTCTRL)

        # Separator radio buttons
        separator_label = wx.StaticText(side_panel, label="Separator:")
        self.separator_choices = [",", ";", ":"]
        self.separator_radio = wx.RadioBox(side_panel, choices=self.separator_choices, majorDimension=1, style=wx.RA_SPECIFY_COLS)

        # Column selection
        col1_label = wx.StaticText(side_panel, label="Select Gene Column")
        self.col1_input = wx.ComboBox(side_panel, choices=[], style=wx.CB_READONLY)
        col2_label = wx.StaticText(side_panel, label="Select P-value Column")
        self.col2_input = wx.ComboBox(side_panel, choices=[], style=wx.CB_READONLY)
        col3_label = wx.StaticText(side_panel, label="Select LogFC Column")
        self.col3_input = wx.ComboBox(side_panel, choices=[], style=wx.CB_READONLY)

        # P-value and LogFC sliders
        pval_logfc_label = wx.StaticText(side_panel, label="Define the p-value and LogFC value to identify DEGs")
        pval_label = wx.StaticText(side_panel, label="P-value:")
        self.pval_slider = wx.Slider(side_panel, value=1, minValue=0, maxValue=100, style=wx.SL_HORIZONTAL)
        log_label = wx.StaticText(side_panel, label="LogFC:")
        self.log_slider = wx.Slider(side_panel, value=1, minValue=0, maxValue=5)
        # Create Static text to display slider value
        self.slider_value_text = wx.StaticText(side_panel, label='0.01', style=wx.ALIGN_CENTER)
        side_panel.SetSizerAndFit(side_sizer)
        self.log_slider_value_text = wx.StaticText(side_panel, label='1', style=wx.ALIGN_CENTER)
        side_panel.SetSizerAndFit(side_sizer)
        # Set default value to 0.01 and 1
        self.pval_slider.SetValue(1)
        self.log_slider.SetValue(1)
        
        self.slider_pval_value = float(self.slider_value_text.GetLabel())
        self.slider_logfc_value = float(self.log_slider_value_text.GetLabel())

        # Expression table
        self.expression_table = wx.grid.Grid(main_panel)
        self.expression_table.CreateGrid(0, 4)
        self.expression_table.SetColLabelValue(0, "Gene")
        self.expression_table.SetColLabelValue(1, "P-value")
        self.expression_table.SetColLabelValue(2, "LogFC")
        self.expression_table.SetColLabelValue(3, "Expression Profile")
        self.expression_table.AutoSizeColumns()
        # Create a save button
        self.save_csv_button = wx.Button(side_panel, label="Save table as CSV")
        
        # Volcano plot
        self.volcano_plot = wx.Panel(main_panel)
        # Create a save button
        self.save_png_button = wx.Button(side_panel, label="Save plot as PNG")
        
        # Clear data button
        self.clear_data_button = wx.Button(side_panel, label="Clear Data")
        self.clear_data_button.SetBackgroundColour("#CD5C5C")
        self.clear_data_button.SetForegroundColour(wx.WHITE)
        
        # Create the "Next" button
        #self.next_button = wx.Button(main_panel, -1, "Next")
        
        # Bind event
        self.file_input.Bind(wx.EVT_FILEPICKER_CHANGED, self.on_upload_button_clicked)
        self.expression_table.Bind(wx.EVT_FILEPICKER_CHANGED, self.update_expression_table)
        self.separator_radio.Bind(wx.EVT_RADIOBOX, self.on_separator_radio_changed)
        self.col1_input.Bind(wx.EVT_COMBOBOX, self.on_col_selected1)
        self.col2_input.Bind(wx.EVT_COMBOBOX, self.on_col_selected2)
        self.col3_input.Bind(wx.EVT_COMBOBOX, self.on_col_selected3)
        self.pval_slider.Bind(wx.EVT_SCROLL, self.on_pval_slider_scroll)
        self.log_slider.Bind(wx.EVT_SCROLL, self.on_logfc_slider_scroll)
        self.disease_combo.Bind(wx.EVT_TEXT, self.on_gene_selection)
        self.save_csv_button.Bind(wx.EVT_BUTTON, self.on_save_csv_button_clicked)
        self.save_png_button.Bind(wx.EVT_BUTTON, self.on_save_png_button_clicked)
        self.clear_data_button.Bind(wx.EVT_BUTTON, self.on_clear_data_button_clicked)
        #self.next_button.Bind(wx.EVT_BUTTON, self.on_next_page)
        #self.col1_input.Bind(wx.EVT_COMBOBOX, self.on_column_selection)

        # Set sizers
        side_panel.SetSizerAndFit(side_sizer)
        main_panel.SetSizerAndFit(main_sizer)
        # Arrange the panels side by side
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(side_panel, 1, wx.EXPAND)
        sizer.Add(main_panel, 5, wx.EXPAND)
        self.SetSizer(sizer)
        
        side_sizer.Add(disease_input_label, 0, wx.ALL | wx.RIGHT, 5)
        side_sizer.Add(self.disease_combo, 0, wx.ALL, 5)
        side_sizer.Add(file_label, 0, wx.ALL, 5)
        side_sizer.Add(self.file_input, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(separator_label, 0, wx.ALL, 5)
        side_sizer.Add(self.separator_radio, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(col1_label, 0, wx.ALL, 5)
        side_sizer.Add(self.col1_input, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(col2_label, 0, wx.ALL, 5)
        side_sizer.Add(self.col2_input, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(col3_label, 0, wx.ALL, 5)
        side_sizer.Add(self.col3_input, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(pval_logfc_label, 0, wx.ALL, 5)
        side_sizer.Add(pval_label, 0, wx.ALL, 5)
        side_sizer.Add(self.pval_slider, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(self.slider_value_text, 0, wx.ALIGN_CENTER)
        side_sizer.Add(log_label, 0, wx.ALL, 5)
        side_sizer.Add(self.log_slider, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(self.log_slider_value_text, 0, wx.ALIGN_CENTER)
        side_sizer.Add(self.save_csv_button, 0, wx.ALIGN_CENTER)
        side_sizer.Add(self.save_png_button, 0, wx.ALIGN_CENTER)
        side_sizer.Add(self.clear_data_button, 0, wx.ALIGN_CENTER)
        
        main_sizer.Add(self.expression_table, 1, wx.EXPAND | wx.ALL)
        #main_sizer.Add(self.next_button, 0, wx.ALIGN_CENTER)
        main_sizer.Add(self.volcano_plot, 2, wx.EXPAND | wx.ALL)
        
    
    def on_gene_selection(self):
        self.biom_list = full_list.loc[full_list['disease'] == self.disease_combo.GetValue(), 'gene']
        
    def on_col_selected1(self,event):
        self.selected_column1 = self.col1_input.GetValue()
        num_rows = self.df.shape[0]
        # Update the column label 
        self.expression_table.SetColLabelValue(0, self.selected_column1)
        # Add the rows to the table based on the selected column
        for i in range(num_rows): 
            self.expression_table.AppendRows()
            self.expression_table.SetCellValue(i, 0, str(self.df.iloc[i][self.selected_column1])) 
        self.expression_table.AutoSizeColumns()
    
    def on_col_selected2(self,event):
        self.selected_column2 = self.col2_input.GetValue()
        num_rows = self.df.shape[0]
        # Update the column label 
        self.expression_table.SetColLabelValue(1, self.selected_column2)
        # Add the rows to the table based on the selected column
        for i in range(num_rows): 
            self.expression_table.AppendRows()
            self.expression_table.SetCellValue(i, 1, str(self.df.iloc[i][self.selected_column2])) 
        self.expression_table.AutoSizeColumns()
    
    def on_col_selected3(self,event):
        self.selected_column3 = self.col3_input.GetValue()
        num_rows = self.df.shape[0]
        # Update the column label 
        self.expression_table.SetColLabelValue(2, self.selected_column3)
        # Add the rows to the table based on the selected column
        for i in range(num_rows): 
            self.expression_table.AppendRows()
            self.expression_table.SetCellValue(i, 2, str(self.df.iloc[i][self.selected_column3])) 
        self.expression_table.AutoSizeColumns()
    
    
    def upload_file(self, file_path):
        # Read data from expression file with the selected separator
        self.df = pd.read_csv(file_path, sep=",", low_memory=False)  ##### the separator is ,
        # Filter genes based on the biomarker list
        self.df = self.df[self.df['gene'].isin(self.biom_list)]
        self.df = self.df.drop_duplicates()
        # Extract desired columns 
        self.df_display = self.df[["gene", "pvalue", "log2FoldChange"]]
        # Update combo boxes
        self.col1_input.Set(self.df.columns)
        self.col2_input.Set(self.df.columns)
        self.col3_input.Set(self.df.columns)
        # Clear existing data in expression table
        self.expression_table.ClearGrid()
        # Update expression table with the extracted column data
        num_rows, num_cols = self.df_display.shape
        self.expression_table.AppendRows(num_rows)
        for row in range(num_rows):
            for col in range(num_cols):
                # Get the cell value 
                cell_value = self.df_display.iloc[row, col]
                # Set the cell value to the corresponding cell in the table
                self.expression_table.SetCellValue(row, col, str(cell_value))

        # Update table layout and create the volcano plot
        self.expression_table.Refresh()
        self.expression_table.AutoSizeColumns()
        self.update_expression_profile()
        self.create_volcano_plot()


    def on_upload_button_clicked(self, event):
        self.on_gene_selection()
        file_path = self.file_input.GetPath()
        self.upload_file(file_path)
        
    def update_expression_table(self,event):
        num_rows, num_cols = self.expression_table.shape  
        # Check if the table is empty
        if self.expression_table.GetNumberRows() == 0:
            self.expression_table.AppendRows(num_rows)
        # Check if the table has a different number of rows
        elif num_rows != self.expression_table.GetNumberRows():
            self.expression_table.ClearGrid()
            
        separator = self.separator_radio.GetStringSelection()

        for row in range(num_rows):
            for col in range(num_cols):
                # Get the cell value
                cell_value = self.df_display.iloc[row, col]
                # Split the cell value using the selected separator
                cell_value_parts = cell_value.split(separator)
                # Set the cell value parts to the corresponding cells in the table
                table_col = col  # Separate variable for table column index
                for i, part in enumerate(cell_value_parts):
                    # Update the table column index to use table_col + i
                    self.expression_table.SetCellValue(row -num_rows, table_col + i, part)

        for row in range(num_rows):
            logfc = float(self.expression_table.GetCellValue(row - num_rows, 2))
            pval = float(self.expression_table.GetCellValue(row - num_rows, 1))
            expression_profile = ""
            if logfc > self.slider_logfc_value and pval < self.slider_pval_value:
                expression_profile = "Upregulated genes"
            elif logfc < -self.slider_logfc_value and pval < self.slider_pval_value:
                expression_profile = "Downregulated genes"
            else:
                expression_profile = "Not Significant"
            self.expression_table.SetCellValue(row - num_rows, 3, expression_profile)

        self.on_gene_selection()
        # Update table layout and volcano plot
        self.expression_table.Refresh()
        self.expression_table.AutoSizeColumns()
        self.update_expression_profile()
        self.create_volcano_plot()
        
        
    def on_separator_radio_changed(self, event):
        separator = self.separator_choices[self.separator_radio.GetSelection()]
        num_rows = self.expression_table.GetNumberRows()
        num_cols = self.expression_table.GetNumberCols()
        data = []
        for row in range(num_rows):
            row_data = []
            for col in range(num_cols):
                row_data.append(self.expression_table.GetCellValue(row, col))
            data.append(row_data)
        # Update expression table with new separator
        df = pd.DataFrame(data, columns=["Gene", "P-value", "LogFC", "Expression Profile"])
        self.update_expression_table(df, separator)
        
        
    def on_pval_slider_scroll(self, event):
        slider = event.GetEventObject()
        self.slider_pval_value = slider.GetValue() / 100  # Convert back to float by dividing by 100
        self.slider_value_text.SetLabel(f'{self.slider_pval_value:.2f}')
        self.update_expression_profile()
    
        
    def on_logfc_slider_scroll(self, event):
        slider = event.GetEventObject()
        self.slider_logfc_value = slider.GetValue()
        self.log_slider_value_text.SetLabel(f'{self.slider_logfc_value:.2f}')
        self.update_expression_profile()
        
    def update_expression_profile(self):
        num_rows = self.expression_table.GetNumberRows()
        for row in range(num_rows):
            logfc = self.expression_table.GetCellValue(row, 2)
            pval = self.expression_table.GetCellValue(row, 1)
            expression_profile = ""
            if self.expression_table.GetColLabelValue(1) == "padj":  
                padj = float(pval)
                if float(logfc) > self.slider_logfc_value and padj < self.slider_pval_value :
                    expression_profile = "Upregulated genes"
                elif float(logfc) < -(self.slider_logfc_value) and padj < self.slider_pval_value :
                    expression_profile = "Downregulated genes"
                else:
                    expression_profile = "Not Significant"
            else:
                if float(logfc) > self.slider_logfc_value and float(pval) < self.slider_pval_value :
                    expression_profile = "Upregulated genes"
                elif float(logfc) < -(self.slider_logfc_value) and float(pval) < self.slider_pval_value :
                    expression_profile = "Downregulated genes"
                else:
                    expression_profile = "Not Significant"
            self.expression_table.SetCellValue(row, 3, expression_profile)
        self.expression_profile_data = [self.expression_table.GetCellValue(row, 3) for row in range(num_rows)]
        self.create_volcano_plot(self.expression_profile_data)
        
    def on_save_csv_button_clicked(self, event):
        # Ask the user where to save the CSV file
        save_dialog = wx.FileDialog(self, "Save as CSV", wildcard="CSV files (*.csv)|*.csv", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_dialog.ShowModal() == wx.ID_CANCEL:
            return
        # Get the file path from the save dialog
        file_path = save_dialog.GetPath()
        # Save the table as a CSV file
        with open(file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Gene', 'P-value', 'LogFC', 'Expression Profile'])
            for row in range(self.expression_table.GetNumberRows()):
                writer.writerow([
                    self.expression_table.GetCellValue(row, 0),
                    self.expression_table.GetCellValue(row, 1),
                    self.expression_table.GetCellValue(row, 2),
                    self.expression_table.GetCellValue(row, 3),])
    
    def create_volcano_plot(self, expression_profile_data):
        # Clear previous plot
        self.volcano_plot.ClearBackground()
        # Retrieve data from selected rows in expression table
        num_rows = self.expression_table.GetNumberRows()
        pval_data = []
        logfc_data = []
        expression_profile_data = []
        for row in range(num_rows):
            pval = float(self.expression_table.GetCellValue(row, 1))
            logfc = float(self.expression_table.GetCellValue(row, 2))
            expression_profile = self.expression_table.GetCellValue(row, 3)
            pval_data.append(pval)
            logfc_data.append(logfc)
            expression_profile_data.append(expression_profile)
        # Create volcano plot
        fig = Figure(figsize=(10, 8), dpi=80)
        canvas = FigureCanvas(self.volcano_plot, -1, fig)
        ax = fig.add_subplot(111)
        # Create a colormap for mapping expression profiles to colors
        color_map = {'Upregulated genes': 'blue', 'Downregulated genes': 'red', 'Not Significant': 'green'}
        colors = np.array([color_map.get(profile, 'black') for profile in expression_profile_data])
        # Plot the data points with colors
        ax.scatter(logfc_data, -np.log10(pval_data), c=colors, edgecolors='none', alpha=1)
        # Add legend
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[profile], markersize=8) for profile in color_map]
        ax.legend(legend_elements, color_map.keys(), title='Expression Profile')
        ax.set_xlabel('LogFC')
        ax.set_ylabel('-log10(P-value)')
        ax.set_title('Volcano Plot')
        # Update the plot
        self.canvas = FigureCanvas(self.volcano_plot, -1, fig)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(canvas, 1, wx.EXPAND)
        self.volcano_plot.SetSizer(self.sizer)
        self.volcano_plot.Layout()
        return fig
        
    def on_save_png_button_clicked(self, event):
        # Ask the user where to save the PNG file
        save_dialog = wx.FileDialog(self, "Save as PNG", wildcard="PNG files (*.png)|*.png", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_dialog.ShowModal() == wx.ID_CANCEL:
            return
        # Get the file path from the dialog
        filepath = save_dialog.GetPath()
        save_dialog.Destroy()
        # Create a new figure for saving as PNG
        fig = self.create_volcano_plot(self.expression_profile_data)
        fig.savefig(filepath, format='png', dpi=300)    
        
    def on_clear_data_button_clicked(self, event):
        # Clear data in expression table
        num_rows = self.expression_table.GetNumberRows()
        self.expression_table.ClearGrid()
        self.expression_table.DeleteRows(0, num_rows)
        self.expression_table.AutoSizeColumns()
        # Clear file input
        self.file_input.SetPath("")
        # Delete the plot
        if self.canvas is not None:
            self.canvas.Destroy()
            self.canvas = None
        self.volcano_plot.DestroyChildren()
        self.volcano_plot.Layout()
    
        self.Layout()
        
#######################################################################################
###################Tab3: Variants######################################################      
#######################################################################################

class TabThree(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        # Create panels
        side_panel = wx.Panel(self)
        side_panel.SetBackgroundColour('#D3D3D3')
        main_panel = wx.Panel(self)
        side_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        #filter by biomarker by disease
        diseases = full_list['disease'].unique()
        disease_input_label = wx.StaticText(side_panel, label='Select the disease of interest:')
        self.disease_combo = wx.ComboBox(side_panel, choices=list(diseases), style=wx.CB_READONLY)
        # File Folder 
        folder_path_label = wx.StaticText(side_panel, label="Enter folder path:")
        self.folder_path_input = wx.TextCtrl(side_panel)
        folder_path_note = wx.StaticText(side_panel, label="Note: Make sure the path contains two folders named 'control' and 'patient'")
        self.variant_type_radio = wx.RadioBox(side_panel, choices=["SNP", "Indels"], majorDimension=1, style=wx.RA_SPECIFY_COLS)
        self.submit_button = wx.Button(side_panel, label="Submit")
        self.save_csv_button = wx.Button(side_panel, label="Save table as CSV")
        # Clear data button
        self.clear_data_button = wx.Button(side_panel, label="Clear Data")
        self.clear_data_button.SetBackgroundColour("#CD5C5C")
        self.clear_data_button.SetForegroundColour(wx.WHITE)
        
        
        # Sizers
        side_sizer.Add(disease_input_label, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(self.disease_combo, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(folder_path_label, 0, wx.ALL, 5)
        side_sizer.Add(self.folder_path_input, 0, wx.EXPAND | wx.ALL, 5)
        side_sizer.Add(folder_path_note, 0, wx.ALL, 5)
        side_sizer.Add(self.variant_type_radio, 0, wx.ALL, 5)
        side_sizer.Add(self.submit_button, 0, wx.ALIGN_CENTER)
        side_panel.SetSizerAndFit(side_sizer)
        side_sizer.Add(self.save_csv_button, 0, wx.ALIGN_CENTER)
        side_sizer.Add(self.clear_data_button, 0, wx.ALIGN_CENTER)
        
        
        # Variant table :SNP
        self.snp_table = wx.grid.Grid(main_panel)
        self.snp_table.CreateGrid(0, 8)
        self.snp_table.SetColLabelValue(0, "Chrom")
        self.snp_table.SetColLabelValue(1, "Gene")
        self.snp_table.SetColLabelValue(2, "SNP position")
        self.snp_table.SetColLabelValue(3, "SNP ID")
        self.snp_table.SetColLabelValue(4, "Reference Genome Allele")
        self.snp_table.SetColLabelValue(5, "Control's Allele")
        self.snp_table.SetColLabelValue(6, "Patient's Allele")
        self.snp_table.SetColLabelValue(7, "Compare")
        self.snp_table.AutoSizeColumns()
        main_sizer.Add(self.snp_table, 1, wx.EXPAND | wx.ALL, 5)
        main_panel.SetSizerAndFit(main_sizer)
        
        # Variant table :Indel
        self.indel_table = wx.grid.Grid(main_panel)
        self.indel_table.CreateGrid(0, 9)
        self.indel_table.SetColLabelValue(0, "Chrom")
        self.indel_table.SetColLabelValue(1, "Gene 1")
        self.indel_table.SetColLabelValue(2, "Gene 2")
        self.indel_table.SetColLabelValue(3, "SNP position")
        self.indel_table.SetColLabelValue(4, "ALT Patient")
        self.indel_table.SetColLabelValue(5, "REF Control")
        self.indel_table.SetColLabelValue(6, "ALT Control")
        self.indel_table.SetColLabelValue(7, "REF Patient")
        self.indel_table.SetColLabelValue(8, "Compare")
        self.indel_table.AutoSizeColumns()
        main_sizer.Add(self.indel_table, 1, wx.EXPAND | wx.ALL, 5)
        main_panel.SetSizerAndFit(main_sizer)
       
        # Bind event
        self.disease_combo.Bind(wx.EVT_TEXT, self.on_gene_selection)
        self.submit_button.Bind(wx.EVT_BUTTON, self.on_submit_button_click)
        self.variant_type_radio.Bind(wx.EVT_RADIOBUTTON, self.on_submit_button_click)
        self.save_csv_button.Bind(wx.EVT_BUTTON, self.on_save_csv_button_clicked)
        self.variant_type_radio.Bind(wx.EVT_RADIOBOX, self.on_variant_type_selection)
        self.clear_data_button.Bind(wx.EVT_BUTTON, self.on_clear_data_button_clicked)
        
        # Set up sizer
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(side_panel, 0, wx.EXPAND)
        sizer.Add(main_panel, 1, wx.EXPAND)
        self.SetSizerAndFit(sizer)
        
        # Hide the Indel table by default
        self.indel_table.Hide()
        
        # Set the initial table display based on the selected radio button
        self.on_variant_type_selection(None)
        
        
    def on_gene_selection(self):
        self.biom_list = full_list.loc[full_list['disease'] == self.disease_combo.GetValue(), 'gene']
        
    def on_variant_type_selection(self, event):
        # Determine which table to show based on the selected radio button
        if self.variant_type_radio.GetStringSelection() == "SNP":
            table_to_show = self.snp_table
            table_to_hide = self.indel_table
        else:
            table_to_show = self.indel_table
            table_to_hide = self.snp_table
        
        # Hide the table that's not selected
        table_to_hide.Hide()
        # Show the table that's selected
        table_to_show.Show()
        # Refresh the layout to reflect the changes
        self.Layout()
        
        
    def read_vcf(self, path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t').rename(columns={'#CHROM': 'CHROM'})  
    
    def read_combined_vcf_files(self,folder_path):
        # Define transition (Ti) and transveersion (Tv)
        ti = ["A>G","G>A","C>T","T>C"]
        tv = ["A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G"]
        # Read data from control 
        folder_name1 = "control"
        folder_path1 = os.path.join(folder_path, folder_name1)
        file_paths_control = glob.glob(os.path.join(folder_path1, '*.vcf'))
        # Read and combine files
        file_contents_control = []
        for path in file_paths_control:
            vcf_reader_control = self.read_vcf(path)
            file_contents_control.append(vcf_reader_control)
        # Concatenate all tables
        self.combined_control = file_contents_control[0]
        for df in file_contents_control[1:]:
            self.combined_control = self.combined_control.append(df, ignore_index=True)
        # Create a new column "nuSub" by concatenating REF and ALT columns
        self.combined_control['nuSub'] = self.combined_control['REF'] + '>' + self.combined_control['ALT']
        # Set values in "TiTv" column based on "nuSub" values
        self.combined_control.loc[self.combined_control['nuSub'].isin(ti), 'TiTv'] = 'Ti'
        self.combined_control.loc[self.combined_control['nuSub'].isin(tv), 'TiTv'] = 'Tv'
        # Rename columns "nuSub" and "TiTv" to "snp_controls" and "TiTv_controls"
        self.combined_control = self.combined_control.rename(columns={'nuSub': 'snp_controls', 'TiTv': 'TiTv_controls'})
        # Read data from patient
        folder_name2 = "patient"
        folder_path2 = os.path.join(folder_path, folder_name2)
        file_paths_patient = glob.glob(os.path.join(folder_path2, '*.vcf'))
        # Read and combine files
        file_contents_patient = []
        for path in file_paths_patient:
            vcf_reader_patient = self.read_vcf(path)
            file_contents_patient.append(vcf_reader_patient)
        # Concatenate all tables
        self.combined_patient = file_contents_patient[0]
        for df in file_contents_patient[1:]:
            self.combined_patient = self.combined_patient.append(df, ignore_index=True)
        # Create a new column "nuSub" by concatenating REF and ALT columns
        self.combined_patient['nuSub'] = self.combined_patient['REF'] + '>' + self.combined_patient['ALT']
        # Set values in "TiTv" column based on "nuSub" values
        self.combined_patient.loc[self.combined_patient['nuSub'].isin(ti), 'TiTv'] = 'Ti'
        self.combined_patient.loc[self.combined_patient['nuSub'].isin(tv), 'TiTv'] = 'Tv'
        # Rename columns "nuSub" and "TiTv" to "snp_controls" and "TiTv_controls"
        self.combined_patient = self.combined_patient.rename(columns={'nuSub': 'snp_patients', 'TiTv': 'TiTv_patients'})
        #print(combined_patient) 
        
    def read_compared_group_final(self,combined_patient,combined_control):
        # Merge patient and control tables based on "CHROM" and "POS" columns
        self.compare_group = pd.merge(combined_patient, combined_control, on=["CHROM", "POS"], how="outer", suffixes=('_patient', '_control'))
        # Replace NaN values with "no"
        self.compare_group = self.compare_group.fillna("no")
        # Set values in "compare" column based on conditions
        self.compare_group.loc[(self.compare_group['snp_controls'] != "no") & (self.compare_group['snp_patients'] == "no"), 'compare'] = "deletion"
        self.compare_group.loc[(self.compare_group['snp_controls'] == "no") & (self.compare_group['snp_patients'] != "no"), 'compare'] = "addition"
        self.compare_group.loc[(self.compare_group['snp_controls'] == self.compare_group['snp_patients']), 'compare'] = "population specific"
        self.compare_group.loc[(self.compare_group['snp_controls'] != self.compare_group['snp_patients']) & (self.compare_group['snp_controls'] != "no") & (self.compare_group['snp_patients'] != "no"), 'compare'] = "different"
        # Free up memory space
        #del file_contents_control
        #del file_contents_patient
        #del combined_control
        #del combined_patient
        
    def read_combined_indels(self,folder_path):
        # Read data from control 
        folder_name1 = "control"
        folder_path1 = os.path.join(folder_path, folder_name1)
        file_paths_control = glob.glob(os.path.join(folder_path1, '*.vcf'))
        # Read and combine files
        file_contents_control = []
        for path in file_paths_control:
            vcf_reader_control = self.read_vcf(path)
            file_contents_control.append(vcf_reader_control)
        # Concatenate all
        self.combined_control = file_contents_control[0]
        for df in file_contents_control[1:]:
            self.combined_control = self.combined_control.append(df, ignore_index=True)
        self.combined_control = self.combined_control.iloc[:, :5]
        self.combined_control = self.combined_control.rename(columns={"#CHROM":"chrom","POS":"POS_snp_control","ID":"id_control","REF":"REF_control","ALT":"ALT_control"})
        # Read data from patient
        folder_name2 = "patient"
        folder_path2 = os.path.join(folder_path, folder_name2)
        file_paths_patient = glob.glob(os.path.join(folder_path2, '*.vcf'))
        # Read and combine files
        file_contents_patient = []
        for path in file_paths_patient:
            vcf_reader_patient = self.read_vcf(path)
            file_contents_patient.append(vcf_reader_patient)
        # Concatenate all tables
        self.combined_patient = file_contents_patient[0]
        for df in file_contents_patient[1:]:
            self.combined_patient = self.combined_patient.append(df, ignore_index=True)
        self.combined_patient = self.combined_patient.iloc[:, :5]
        self.combined_patient = self.combined_patient.rename(columns={"#CHROM":"chrom","POS":"POS_snp_patient","ID":"id_patient","REF":"REF_patient","ALT":"ALT_patient"})
            
            
    def on_save_csv_button_clicked(self, event):
        # Ask the user which table to save
        table_choice = wx.SingleChoiceDialog(self, "Which data do you want to save?", "Choose table", ["Indels", "SNPs"])
        if table_choice.ShowModal() == wx.ID_CANCEL:
            return
        table_name = table_choice.GetStringSelection()

        # Ask the user where to save the CSV file
        save_dialog = wx.FileDialog(self, "Save as CSV", wildcard="CSV files (*.csv)|*.csv", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_dialog.ShowModal() == wx.ID_CANCEL:
            return 

        # Get the file path from the save dialog
        file_path = save_dialog.GetPath()
        # Define the column names for both tables
        indel_cols = ['Chrom', 'Gene1', 'Gene2', ' SNP position', 'ALT Patient',"REF Control","ALT Control","REF Patient" 'Compare']
        snp_cols = ['Chrom', 'Gene', 'SNP position', 'SNP ID', 'Reference Genome Allele', 'Control\'s Allele', 'Patient\'s Allele', 'Compare']

        # Save the selected table data as a CSV file
        with open(file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            if table_name == "Indels":
                writer.writerow(indel_cols)
                table = self.indel_table
            else:
                writer.writerow(snp_cols)
                table = self.snp_table
            for row in range(table.GetNumberRows()):
                writer.writerow([
                    table.GetCellValue(row, 0),
                    table.GetCellValue(row, 1),
                    table.GetCellValue(row, 2),
                    table.GetCellValue(row, 3),
                    table.GetCellValue(row, 4),
                    table.GetCellValue(row, 5),
                    table.GetCellValue(row, 6),
                    table.GetCellValue(row, 7)]) 
        

    def on_submit_button_click(self, event):
        # display the message dialog
        msg = wx.BusyInfo("Processing data, please wait...")
        # Get input values
        folder_path = self.folder_path_input.GetValue()
        vtype = self.variant_type_radio.GetSelection()
        if vtype == 0 :
            # Read VCF as df
            self.read_combined_vcf_files(folder_path)
            self.read_compared_group_final(self.combined_patient,self.combined_control)
            # Get data from compare_group
            self.compare_group_final = self.compare_group.iloc[:, [0, 1, 2, 3, 4, 13, 14, 17, 25, 26, 27]]
            self.compare_group_final.columns = ["Chrom"] + self.compare_group_final.columns.tolist()[1:]
            self.compare_group_final = self.compare_group_final.replace(to_replace="chr", value="", regex=True)
            # Extract columns Chrom and POS from compare_group_final
            chrom_snppos = self.compare_group_final.iloc[:, [0, 1]]
            chrom_snppos.columns = ["chrom", "snp_position"]
            chrom_snppos = chrom_snppos.replace(to_replace="chr", value="", regex=True)
            # Extract gene positions
            gene_pos = annotation1.iloc[:, [4, 5, 6, 8]].drop_duplicates()
            gene_pos.columns = ["chrom", "start", "end", "gene name"]
            gene_pos['start'] = gene_pos['start'].astype('int32')
            gene_pos['end'] = gene_pos['end'].astype('int32')
            # Merge gene_pos and chrom_snppos
            ann_snps = pd.merge(gene_pos, chrom_snppos, on="chrom")
            del gene_pos
            del chrom_snppos
            ann_snps = ann_snps[(ann_snps["snp_position"] < ann_snps["end"]) & (ann_snps["snp_position"] > ann_snps["start"])].reset_index(drop=True)
            ann_snps = ann_snps.drop(columns=["end"])
            # Merge ann_snps and compare_group_final
            ann_snps2 = pd.merge(ann_snps, self.compare_group_final, left_on=["chrom", "snp_position"], right_on=["Chrom", "POS"])
            del ann_snps
            del self.compare_group_final
            ann_snps3 = ann_snps2.iloc[:, [0, 2, 3, 6, 7, 11, 8, 14]]
            ann_snps3.columns = ["Chrom", "Gene", "SNP position", "SNP ID", "Reference Genome Allele", "Control's Allele", "Patient's Allele", "Compare"]
            ann_snps3 = ann_snps3.replace(to_replace="no", value=".", regex=True)
            ann_snps3 =ann_snps3.drop_duplicates()
            #print("ann_snps3:")
            #print(ann_snps3)
            self.on_gene_selection()
            final_table = ann_snps3[ann_snps3['Gene'].isin(self.biom_list)]
            #print(final_table)
            # Update SNP table with ann_snps3 table data
            self.snp_table.ClearGrid()
            num_rows, num_cols = final_table.shape
            self.snp_table.AppendRows(num_rows)
            for row in range(num_rows):
                for col in range(num_cols):
                    # Get the cell value 
                    cell_value = final_table.iloc[row, col]
                    # Set the cell value to the corresponding cell in the table
                    self.snp_table.SetCellValue(row, col, str(cell_value))
            # Update table layout
            self.snp_table.Refresh()
            self.snp_table.AutoSizeColumns()
            
        elif vtype == 1:
            # Read VCF as df
            self.read_combined_indels(folder_path)
            # Extract columns chrom and snp position
            chrom_snppos = self.combined_control.iloc[:, [0, 1]]
            chrom_snppos.columns = ["chrom", "snp_position"]
            chrom_snppos = chrom_snppos.replace(to_replace="chr", value="", regex=True)
            # Extract gene positions
            gene_pos = annotation1.iloc[:, [4, 5, 6, 8]].drop_duplicates()
            gene_pos.columns = ["chrom", "start", "end","gene name"]
            # Merge gene_pos and chrom_snppos
            annotated_indel = pd.merge(gene_pos, chrom_snppos, on="chrom")
            del gene_pos
            del chrom_snppos
            annotated_indel = annotated_indel[(annotated_indel["snp_position"] < annotated_indel["end"]) & (annotated_indel["snp_position"] > annotated_indel["start"])].reset_index(drop=True)
            annotated_indel = annotated_indel.drop(columns=["end"])
            # Remove 'chr' from chromosome column in combined_control/patient tables
            self.combined_control = self.combined_control.replace(to_replace="chr", value="", regex=True)
            self.combined_patient = self.combined_patient.replace(to_replace="chr", value="", regex=True)
            #print("self.combined_control")
            #print(self.combined_control)
            #print("self.combined_patient")
            #print(self.combined_patient)
            #print("annotated_indel")
            #print(annotated_indel)
            # Merge annotated_indel with combined_control on matching chrom and snp_position
            annotated_indel2 = pd.merge(annotated_indel, self.combined_control, left_on=['chrom', 'snp_position'], right_on=['CHROM', 'POS_snp_control'])
            annotated_indel2 = annotated_indel2.drop(columns=["CHROM" , "POS_snp_control"])
            #print("annotated_indel2")
            #print(annotated_indel2)

            # Extract chrom and snp_position columns from combined_patients_indel and remove 'chr' from chrom values
            self.combined_patient.rename(columns={"CHROM": "chrom"}, inplace=True)
            chrom_snppos = self.combined_patient[['chrom', 'POS_snp_patient']]
            chrom_snppos = chrom_snppos.replace(to_replace="chr", value="", regex=True)
            #print("chrom_snppos")
            #print(chrom_snppos)

            # Extract chrom, start, end, and gene columns from annotation1 and rename columns
            gene_pos = annotation1.iloc[:, [4, 5, 6, 8]].drop_duplicates()
            gene_pos.columns = ["chrom", "start", "end", "gene name"]
            #print("gene_pos")
            #print(gene_pos)

            # Merge gene_pos with chrom_snppos on matching chrom and snp_position with snp_position between gene_start and gene_end
            annotated_indel3 = pd.merge(gene_pos, chrom_snppos, on="chrom")
            #print("annotated_indel3")
            #print(annotated_indel3)
            annotated_indel3 = annotated_indel3[(annotated_indel3['POS_snp_patient'] > annotated_indel3['start']) & (annotated_indel3['POS_snp_patient'] < annotated_indel3['end'])]
            annotated_indel3 = annotated_indel3.drop('POS_snp_patient', axis=1)
            # Merge annotated_indel with combined_patients_indel on matching chrom and snp_position
            annotated_indel4 = pd.merge(annotated_indel, self.combined_patient, left_on=['chrom', 'snp_position'], right_on=['chrom', 'POS_snp_patient'])
            annotated_indel4 = annotated_indel4.drop(['POS_snp_patient'], axis=1)
            annotated_indel4 = annotated_indel4.rename(columns={'gene name': 'gene1'})
            #print("annotated_indel4")
            #print(annotated_indel4)

            # Rename columns in annotated_indel2
            annotated_indel2 = annotated_indel2.rename(columns={'gene name': 'gene2'})
            # Compare groups and assign 'compare' values
            compare_group_i = pd.merge(annotated_indel4, annotated_indel2, how='outer', on=['chrom', 'snp_position'])
            compare_group_i['compare'] = ''
            #print("compare_group_i")
            #print(compare_group_i)
            compare_group_i.loc[(compare_group_i['ALT_control'].notna()) & (compare_group_i['ALT_patient'].isna()), 'compare'] = 'deletion'
            compare_group_i.loc[(compare_group_i['ALT_control'].isna()) & (compare_group_i['ALT_patient'].notna()), 'compare'] = 'addition'
            compare_group_i.loc[(compare_group_i['ALT_control'] == compare_group_i['ALT_patient']) & (compare_group_i['ALT_control'].notna()) & (compare_group_i['ALT_patient'].notna()), 'compare'] = 'population specific'
            compare_group_i.loc[(compare_group_i["ALT_control"] != compare_group_i['ALT_patient']) & (compare_group_i['ALT_control'] !="no") & (compare_group_i['ALT_patient'] !="no"), 'compare'] = "different"
            compare_group_i = compare_group_i[["chrom", "gene1", "gene2", "snp_position", "ALT_patient", "REF_control", "ALT_control", "REF_patient","compare"]]
            compare_group_i.fillna(".", inplace=True)
            compare_group_i =compare_group_i.drop_duplicates()
            self.on_gene_selection()
            final_table = pd.merge(self.biom_list, compare_group_i, how='inner', left_on='gene', right_on='gene1')
            final_table = pd.concat([final_table, pd.merge(self.biom_list, compare_group_i, how='inner', left_on='gene', right_on='gene2')])
            final_table = final_table.drop(['gene'], axis=1)
            #print(final_table)

            # Update Indel table with ann_snps3 table data
            self.indel_table.ClearGrid()
            num_rows, num_cols = final_table.shape
            self.indel_table.AppendRows(num_rows)
            for row in range(num_rows):
                for col in range(num_cols):
                    # Get the cell value
                    cell_value = final_table.iloc[row, col]
                    # Set the cell value to the corresponding cell in the table
                    self.indel_table.SetCellValue(row, col, str(cell_value))
            # Update table layout
            self.indel_table.Refresh()
            self.indel_table.AutoSizeColumns()
            
    def on_clear_data_button_clicked(self, event):
        # Clear data in expression table
        vtype = self.variant_type_radio.GetSelection()
        if vtype == 0 :
            num_rows = self.snp_table.GetNumberRows()
            self.snp_table.ClearGrid()
            self.snp_table.DeleteRows(0, num_rows)
            self.snp_table.AutoSizeColumns()
        elif vtype == 1 :
            num_rows = self.indel_table.GetNumberRows()
            self.indel_table.ClearGrid()
            self.indel_table.DeleteRows(0, num_rows)
            self.indel_table.AutoSizeColumns()
        # Clear file input
        self.folder_path_input.SetValue("")
            
        self.Layout()

class NeuroVar(wx.Frame):
    def __init__(self):
        super().__init__(parent=None, title='NeuroVar')
        # Create a panel and notebook (tabs holder)
        p = wx.Panel(self)
        nb = wx.Notebook(p)
        # Create the tab windows
        Biomarker = TabOne(nb)
        Expression = TabTwo(nb)
        Variants = TabThree(nb)
        # Set the font size for the navigation bar tabs
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        nb.SetFont(font)
        # Add the windows to tabs and name them.
        nb.AddPage(Biomarker, "Biomarker")
        nb.AddPage(Expression, "Expression")
        nb.AddPage(Variants, "Variants")
        # Set noteboook in a sizer to create the layout
        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)
        p.SetSizer(sizer)
        

if __name__ == '__main__':
    app = wx.App()
    frame = NeuroVar()
    #frame.Center()
    screen_width, screen_height = wx.GetDisplaySize()
    frame.SetSize(screen_width, screen_height)
    frame.SetPosition((0, 0))
    frame.Show()
    app.MainLoop()
