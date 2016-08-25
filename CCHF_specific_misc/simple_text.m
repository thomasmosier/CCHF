function simple_text(path,text)

fid = fopen(path,'a');

fprintf(fid,'%s',text);

fclose(fid);

