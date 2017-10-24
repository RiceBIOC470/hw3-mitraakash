% GB Comments
1.	100
2a. 70 Be careful with the use of the function showalignment. Feeding your align into the function only gives a snippet of the entire possible coding sequence and therefore outputs a artificially high percent alignment. If you determine the length of ERK1.sequence, you will notice that its length is not equal to the length of the sequence used in alignment. 
2b. 70 same issue as 2a
2c. 70 same issue as 2a. 
3a 100 
3b. 100
3c. 100  	
Overall: 87


% Akash Mitra
% am132

%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

seq1 = 'GTAATCC';
seq2 = 'GTATCCG';
Score_matrix = [2 -1 -1 -1; -1 2 -1 -1; -1 -1 2 -1; -1 -1 -1 2];
[score1, align1, start1] = swalign(seq1,seq2, 'Alphabet', 'nt' , 'ScoringMatrix' , Score_matrix, 'GapOpen', 1, 'ShowScore', true);

disp(["Score of highest scoring alignment : " int2str(score1)]); 
disp(align1);

% Highest scoring alignment ignores gap at the last position, resulting in a score of 11.
% If you include gap at the end of the alignment, highest score is 10.

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

erk1_info = getgenbank('NM_002746');
erk2_info = getgenbank('NM_002745');

erk1_location = erk1_info.CDS.indices;
erk1_coding_dna = erk1_info.Sequence(erk1_location(1) : erk1_location(2));
erk2_location = erk2_info.CDS.indices;
erk2_coding_dna = erk2_info.Sequence(erk2_location(1) : erk2_location(2));

[score2, align2, start2] = swalign(erk1_coding_dna,erk2_coding_dna, 'Alphabet', 'nt' , 'ShowScore', true);
showalignment(align2);

% 811/1073 = 76% identity alignment b/w two sequences

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

erk1_protein = erk1_info.CDS.translation;
erk2_protein = erk2_info.CDS.translation;

[score3, align3, start3] = swalign(erk1_protein,erk2_protein, 'Alphabet', 'aa' , 'ShowScore', true);
showalignment(align3);

% 305/346 = 88% identity b/w two sequences and 96% positive (similar aa)

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

% Working with mRNA and coding DNA sequences

mouse_erk1_info = getgenbank('NM_011952');
mouse_erk2_info = getgenbank('NM_001038663');

mouse_erk1_location = mouse_erk1_info.CDS.indices;
mouse_erk1_coding_dna = mouse_erk1_info.Sequence(mouse_erk1_location(1) : mouse_erk1_location(2));
mouse_erk2_location = mouse_erk2_info.CDS.indices;
mouse_erk2_coding_dna = mouse_erk2_info.Sequence(mouse_erk2_location(1) : mouse_erk2_location(2));

[score4, align4, start4] = swalign(mouse_erk1_coding_dna,erk1_coding_dna, 'Alphabet', 'nt' , 'ShowScore', true);
showalignment(align4);

% 1026/1137 = 90% identity for coding regions in mouse mRNA and human coding DNA sequence of
% ERK1

[score5, align5, start5] = swalign(mouse_erk2_coding_dna,erk2_coding_dna, 'Alphabet', 'nt' , 'ShowScore', true);
showalignment(align5);

% 996/1075 = 93% identity for coding regions in mouse mRNA and human coding DNA sequence of
% ERK2

% Working with mRNA and protein sequences

mouse_erk1_protein = mouse_erk1_info.CDS.translation;
mouse_erk2_protein = mouse_erk2_info.CDS.translation;

[score6, align6, start6] = swalign(mouse_erk1_protein,erk1_protein, 'Alphabet', 'aa' , 'ShowScore', true);
showalignment(align6);

% 97% identity and 98% positives (similar aa) for mouse translated protein and human protein sequence of
% ERK1

[score7, align7, start7] = swalign(mouse_erk2_protein,erk2_protein, 'Alphabet', 'aa' , 'ShowScore', true);
showalignment(align7);

% 99% identity and 100% positives (similar aa) for mouse translated protein and human protein sequence of
% ERK2

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits.

ans_hit = genbankhits('AF308602', 23);

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

sourceGenbank('XM_001172077');

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

human_gene = 'NM_012271'; % human SETD2
dog_gene = 'NM_001003181'; % dog KIT

sourceGenbank(human_gene);

%Closest match results in the same gene SETD2 hit from humans
%Closest non human match results in a chimpanzee hit

sourceGenbank(dog_gene);

%Human hit not found.
%Closest non human match results in a dog KIT gene hit.
