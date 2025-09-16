%command_output = ssh2_simple_command('128.141.221.114','GDD','Win_Admin','pscp -pw  -r c:\Users\GDD\Picosec fbrunbau@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Pool1')
%command_output = ssh2_simple_command('128.141.221.91','GDD','Win_Admin','pscp -pw  -r d:\Picosec fbrunbau@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Pool1')
javaaddpath('ganymed-ssh2-build250.jar')


ssh2_conn = ssh2_config('128.141.41.251','GDD','Win_Admin')
%ssh2_conn = ssh2_config('128.141.221.91','GDD','Win_Admin')
ssh2_conn = ssh2_command(ssh2_conn, 'pscp -r -pw Win_Admin D:\PicosecTestBeam\Run241 gdd@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Pool1')
%ssh2_conn = ssh2_close(ssh2_conn); %close connection when done


% channel  =  sshfrommatlab('GDD','128.141.221.114','Win_Admin')

% [channel, result]  =  sshfrommatlabissue(channel,'pscp -pw  -r c:\Users\GDD\Picosec fbrunbau@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Pool1')

 
%  channel  =  sshfrommatlab('GDD','128.141.41.251','Win_Admin')

 %[channel, result]  =  sshfrommatlabissue(channel,'pscp -r -pw Win_Admin D:\PicosecTestBeam\Run241 gdd@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Pool1')
