

fid = fopen('log2.html','a');
fprintf(fid, '<H2>Inputs</H2>');
fprintf(fid, ['\n<BR><U>Input String:</U> ', ...
'Stephan kaas']);
fprintf(fid, ['\n<BR><U>Input Boolean:</U> ', ...
num2str(3)]);
fprintf(fid, ['\n<BR><U>Input Number:</U>   '  num2str(4)]);

fprintf(fid, '\n<TABLE><TR><TD><img src="fig2.jpg" alt="no picture" title="SSTUUSS" width=200 height=200></TD></TR></TABLE>' );
fclose(fid);

