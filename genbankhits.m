function x = genbankhits(accession, hits)

if ~exist('hits', 'var')
     hits = 50;
end

seq_info = getgenbank(accession);
[blast_ID, blast_time] = blastncbi(seq_info, 'blastn'); 
blast_report = getblast(blast_ID, 'WaitTime', blast_time);
accessions_name = {};

for i=1:hits
    temp_name = blast_report.Hits(i).Name;
    location_name = strfind(temp_name, '|');
    accession_name = temp_name(location_name(3)+1:location_name(4)-1);
    accessions_name{i} = accession_name;
end

x = accessions_name;

end


