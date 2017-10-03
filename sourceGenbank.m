function sourceGenbank(accession)

accession_hits = genbankhits(accession);

human_hit = [];
nonhuman_hit = [];
nonhuman_org = [];

for i=1:length(accession_hits)
    seq_info = getgenbank(char(accession_hits(i)));   
    if (strcmp(seq_info.Source,'Homo sapiens (human)'))
        human_hit = seq_info.Accession;
        break;
    end
end

for i=1:length(accession_hits)
    seq_info = getgenbank(char(accession_hits(i))); 
    if ~(strcmp(seq_info.Source,'Homo sapiens (human)'))
        nonhuman_hit = seq_info.Accession;
        nonhuman_org = seq_info.Source;
        break;
    end 
end

if ~isempty(human_hit)
    disp(['Closest human match is: ' human_hit]);
else
    disp('No human match found'); 
end

disp(['Closest non human match is: ' nonhuman_hit]);
disp(['Organism is: ' nonhuman_org]);