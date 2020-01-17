mtBreak = read.csv('../../Body/1Raw/MitoBreakDB_12122019.csv')

konstantinData1 = read.table('../../Body/1Raw/compare.square.200.38windows.txt', header=TRUE, sep='\t')
konstantinData2 = read.table('../../Body/1Raw/compare.square.200.38windows_2.txt', header=TRUE, sep='\t')

row.names(konstantinData1) = konstantinData1[, 1]
konstantinData1 = konstantinData1[, -1]

row.names(konstantinData2) = konstantinData2[, 1]
konstantinData2 = konstantinData2[, -1]

a = konstantinData2 / konstantinData1
