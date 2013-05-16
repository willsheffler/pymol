
df  = read_csv(filename, index_col = 0,header = 0)
self.datatable = QtGui.QTableWidget(parent=self)
self.datatable.setColumnCount(len(df.columns))
self.datatable.setRowCount(len(df.index))
for i in range(len(df.index)):
    for j in range(len(df.columns)):
        self.datatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iget_value(i, j))))