classdef suitesparse_helper
    
    methods (Static = true)

        function ids = get(criteria, names)
            % Get index of the SuiteSparse Matrix Collection
            addpath('SuiteSparse-7.1.0/ssget');
            index = ssget;

            % filter based on `criteria`
            size_conditions = (criteria.n_max >= index.ncols) & (index.ncols >= criteria.n_min);
            is_symmetric = index.numerical_symmetry == criteria.symmetric & index.pattern_symmetry == criteria.symmetric;
            conditions = size_conditions & is_symmetric;
            conditions = conditions & index.posdef == criteria.posdef & index.isReal == criteria.real;
            
            n_names = length(names);
            if n_names > 0
                include_names = zeros(size(conditions));
                % filter based on names`
                for i = 1:n_names
                    include_names = include_names | strcmp(index.Name, names(i))';  % (transpose is important!) 
                end
                conditions = conditions & include_names;
            end
            
            % find and return matrices satisfying all conditions
            ids = find(conditions);

            % sort by size
            [~, i] = sort (index.nnz(ids)) ;
            ids = ids(i);

            fprintf('Retrieved %s matrices from SuiteSparse.\n', num2str(length(ids)));
        end
    end
end