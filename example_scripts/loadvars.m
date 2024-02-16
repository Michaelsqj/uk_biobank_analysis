% Function to load variables from a .mat file into the calling function workspace
function loadvars(varnames,fname)

% Force input to be a cell array, splitting by rows if a character array is
% provided
if ~iscell(varnames)
    varnames=mat2cell(varnames,ones(size(varnames,1),1));
end

% Loop through the variables to load
for ii = 1:length(varnames)
    
    % Check if it already exists in the calling workspace
    if evalin('caller', ['exist(''' varnames{ii} ''',''var'')'])
        disp(['Variable ' varnames{ii} ' already exists in the calling workspace, skipping'])
    else
        disp(['Loading variable: ',varnames{ii}])
        
        % Load into a temporary variable
        tmp = load(fname,varnames{ii});
        
        % Assign back to the calling workspace
        % NB. load above creates a struct, so var1 would be loaded as tmp.var1
        % Use eval to undo this
        disp(['Assigning ',varnames{ii} ' to calling workspace'])
        assignin('caller',varnames{ii},eval(['tmp.' varnames{ii}]));
    end
end