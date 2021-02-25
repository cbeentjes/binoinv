% Quick helper function which converts many outputs from
% NormalApproximation.nb into a single Matlab struct
names = {'Error','M','p','RelativeError','W','xs','delta'};
flag = true;
for ix = 1:length(names)
    if isfile(['mathematicadata',names{ix},'.mat'])
        mathematicadata.(names{ix}) = load(['mathematicadata',names{ix},'.mat']);
        mathematicadata.(names{ix}) = mathematicadata.(names{ix}).Expression1;
        delete(['mathematicadata',names{ix},'.mat']);
    else
        flag = false;
    end
end
if flag
    save('mathematicadata.mat','mathematicadata')
end
