function subject_list = pl_experiment_OASIS_subjects(subject_file)
% PL_EXPERIMENT_OASIS_SUBJECTS compiles the filenames of all OASIS subjects
% for which we have segmentations.
%
% Author(s): Roland Kwitt, 2015

subjects = importdata(subject_file);
subject_list =  cell(length(subjects),1);    
    
for subject_id=1:length(subjects)
    parts = strsplit(subjects{subject_id},'_');
    subject_list{subject_id} = sprintf('%s.long.%s_%s_CC.stl', ...
        subjects{subject_id}, parts{1}, parts{2});
end