classdef SuitesSparseHelper
    
    methods (Static = true)

        function ids = get(criteria)
            % Get index of the SuiteSparse Matrix Collection
            addpath('SuiteSparse/ssget');
            index = ssget;

            % filter based on `criteria`
            if isfield(criteria, 'names')
                % either we provide explicit names of desired matrices...
                n_names = length(criteria.names);
                conditions = zeros(size(index.nentries));
                % filter based on `names`
                for i = 1:n_names
                    conditions = conditions | strcmp(index.Name, criteria.names(i))';  % (transpose is important!) 
                end
            else
                % ... or more general criteria
                if ~isfield(criteria, 'n_max')
                    criteria.n_max = inf;
                end
                if ~isfield(criteria, 'n_min')
                    criteria.n_min = 1;
                end
                conditions = (criteria.n_max >= index.ncols) & (index.ncols >= criteria.n_min);
                if isfield(criteria, 'symmetric')
                    is_symmetric = index.numerical_symmetry == criteria.symmetric & index.pattern_symmetry == criteria.symmetric;
                    conditions = conditions & is_symmetric;
                end
                if isfield(criteria, 'posdef')
                    conditions = conditions & index.posdef == criteria.posdef;
                end
                if isfield(criteria, 'real')
                    conditions = conditions & index.isReal == criteria.real;
                end
            end
            if isfield(criteria, 'skip')
                n_skip = length(criteria.skip);
                for i = 1:n_skip
                    conditions = conditions & ~strcmp(index.Name, criteria.skip(i))';  % (transpose is important!) 
                end
            end

            % find and return matrices satisfying all conditions
            ids = find(conditions);

            % sort by size
            [~, i] = sort (index.nnz(ids)) ;
            ids = ids(i);

            fprintf('Found %s matrices.\n', num2str(length(ids)));
        end
    end
end