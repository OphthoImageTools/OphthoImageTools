function batchCorneaProcessingAndAnalysis(parentdir)

list = dir([parentdir filesep '**/*.*']);
subdir = list(cell2mat({list.isdir}));
useful = ~contains({subdir.name},'.');
subdir = subdir(useful);
resultsfile = fullfile(parentdir,['results' date '.csv']);

n = 0;

for i = 1:numel(subdir)
    folder = fullfile(subdir(i).folder, subdir(i).name);
    xml_files = dir(fullfile(folder,'*.xml'));
    precipitatesVolumefile = fullfile(folder,'precipitatesVolume.mat');
    alreadyProcessed = exist(precipitatesVolumefile,'file');

    % if has already been preocessed (i.e. results file in the folder) take
    % the results from the file and add them to the list
    if alreadyProcessed

        n = n+1;
        disp(['Taking results found in folder ' subdir(i).folder filesep subdir(i).name]);

        % generate an ID from the path
        currentParent = subdir(i).folder;
        if length(subdir(i).folder) > length(parentdir)
            ID = [currentParent(length(parentdir)+1:end) '_' subdir(i).name] ;
        else
            ID = [subdir(i).name] ;
        end

        % take saved results from the precipitatesVolumefile
        load(precipitatesVolumefile)
        results(n).ID = ID;
        results(n).precipitates_volume_mm3 = precVol;
        results(n).area_mm2 = area;
        results(n).ratio_PrecVol_area_mm = ratio_PrecVol_area;

    % if there is a single xml file and the folder hansn't already been
    % preocessed (i.e. results file are not in the folder)
    elseif numel(xml_files) == 1 && ~alreadyProcessed

        disp(['Working on xml in folder ' subdir(i).folder filesep subdir(i).name]);

        % process
        n = n+1;
        [precVol, area, ratio_PrecVol_area] = corneaProcessingAndAnalysis(fullfile(subdir(i).folder, subdir(i).name));

        % store results

        % generate an ID from the path
        currentParent = subdir(i).folder;
        if length(subdir(i).folder) > length(parentdir)
            ID = [currentParent(length(parentdir)+2:end) '_' subdir(i).name] ;
        else
            ID = [subdir(i).name] ;
        end
        ID = strrep(ID,'\','_');
        ID = strrep(ID,'/','_');
        ID = strrep(ID,', ','_');
        ID = strrep(ID,'(','_');
        ID = strrep(ID,')','_');
        ID = strrep(ID,' ','_');
        ID = strrep(ID,'__','_');

        results(n).ID = ID;
        results(n).precipitates_volume_mm3 = precVol;
        results(n).area_mm2 = area;
        results(n).ratio_PrecVol_area_mm = ratio_PrecVol_area;

        % save image with ID
        f = getframe(gcf);
        imwrite( f.cdata, fullfile(parentdir, [ID '.png']) );

    end
end

% save results
writetable( struct2table(results), resultsfile);



