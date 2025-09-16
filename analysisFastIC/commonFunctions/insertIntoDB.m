username = "admin";
password = "Enter123!";

vendor = "MySQL";
opts = databaseConnectionOptions("native",vendor);
opts = setoptions(opts, ...
    'DataSourceName',"MySQLDataSource", ...
    'DatabaseName',"assetQE",'Server',"dbod.gc022.cern.ch", ...
    'PortNumber',5512);

status = testConnection(opts,username,password);
saveAsDataSource(opts)

datasource = "MySQLDataSource";
username = "root";
password = "matlab";
conn = mysql(datasource,username,password);

sqlquery = strcat("INSERT INTO picosecAnalysis (dut) VALUES (1)",1);
execute(conn,sqlquery)