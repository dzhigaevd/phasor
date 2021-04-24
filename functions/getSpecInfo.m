function [motname, motval, scan, val_moved, moved_motnames] = getSpecInfo(filename, scanN)
    %SPECN_PARAM Summary of this function goes here
    %   Extract motors and their values from spec file
    fid=fopen(filename);

    F=fread(fid);
    FS=char(F');

    fposS=strfind(FS,'#S');
    scan=[];
    fposS=[fposS length(FS)];
    
    for k=1:length(fposS)-1
        %---- Get first column headers ----
        fseek(fid,fposS(k)-1,'bof');
        com{k}=fgetl(fid);	
        C = strsplit(com{k}(3:8), ' ');
        scan=[scan str2num(C{2})];
        frewind(fid)
    end

    fclose(fid);

    %db=unique(scan);
    %for i=1:length(db)
    %    fdb=find(scan==db(i));
    %    for j=1:length(fdb),scan(fdb(j))=scan(fdb(j))+j/10;end %add .1, .2, ...
    %end

    scan_open = scan(scanN);

    fid=fopen(filename);
    F=fread(fid);FS=char(F');
    fposS1=strfind(FS,['#S 1']);
    fposSc=strfind(FS,['#S ' num2str(floor(scan_open)) ' ']);
    %ranksc=round((scan_open-floor(scan_open))*10);
    %fposSc=fposSc(ranksc);

    %---get motor names-----------------------------
    fposO=strfind(FS,'#O0');
    motname=[];
    fseek(fid,fposO(1)-1,-1);
    while ftell(fid)<min([fposS1])-1
        tmpmot=fgetl(fid);
        motname=[motname tmpmot(4:end) ' '];
    end
    motname=extractstr(motname,'  ');

    %---get motor values-----------------------------
    fseek(fid,fposSc-1,-1);
    tmppos=ftell(fid);
    fposN=strfind(FS(tmppos:end),'#N')+tmppos-1;

    motval=[];
    %fseek(fid,fposP0-1,-1);
    counter=0;
    while ftell(fid)<fposN-1
        tmpmot=fgetl(fid);
        if strmatch('#P',tmpmot) 
            counter=counter+1;
            if counter<10
                motval=[motval tmpmot(4:end) ' '];
            else 
                motval=[motval tmpmot(5:end) ' '];
            end
        end
    end

    motval=extractstr(motval,' ');
    tmpval=[];
    for ii=1:length(motval)
        tmpval=[tmpval str2num(motval{ii})];
    end
    motval=tmpval;

    fseek(fid,fposSc-1,-1);
    tmppos=ftell(fid);
    fposL=strfind(FS(tmppos:end),'#L')+tmppos-1;

    fseek(fid,fposSc+3,-1);
    tmppos=ftell(fid);
    fposSnext=strfind(FS(tmppos:end),'#S')+tmppos-1;

    %%% get moved motors names
    fseek(fid,fposL(1)-1,-1);
    tmppos=ftell(fid);
    moved_motnames_str = fgetl(fid);
    moved_motnames_str = moved_motnames_str(4:end);
    moved_motnames = strsplit(moved_motnames_str,'  ');

    val_moved = [];
    while ftell(fid)<fposSnext-2
        tmpval_moved_str=fgetl(fid);
        val_moved_str = strsplit(tmpval_moved_str,' ');

        if  strcmp(val_moved_str{1},'#C')==0

        for mm = 1:length(val_moved_str)
            val_moved_1(mm) = str2num(val_moved_str{mm});
        end
            val_moved = [val_moved; val_moved_1];
        else
            break
        end
    end

    %%%% get all values for moved motors
    function word=extractstr(s,separator)
        skip=length(separator);
        count=0;

        lf=strfind(s,separator);
        lf=[1-skip lf length(s)+1];

        for jj=1:length(lf)-1
           tmp=s(lf(jj)+skip:lf(jj+1)-1);
           if ~isempty(tmp)
               count=count+1;
               word{count}=tmp;
           end
        end
    end

    specN = [];

    specN.delta = motval(1);
    specN.gamma = motval(6);
    specN.dist = motval(27);
    specN.energy = motval(30);
end