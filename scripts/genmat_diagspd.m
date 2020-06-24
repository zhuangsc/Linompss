## usage: A = genmat_diagspd(m, b, f, inc=0.0) 
## A is m x m, divided into b x b blocks, has a fill approximate to f and inc controls the difference from 1-norm equal to 1.
function A = genmat_diagspd(m, b, f, inc=0.0)

A = sprand(m, m, f);

for i = 1:b:m
	iend = i + b - 1;
	
	# make the diagonal block symmetric and diagonally dominant
	for j=i:iend
		for bi=j+1:iend
			A(bi,j) = A(j,bi);
		end
		A(j,j) += b;
	end

end

for j = 1:m
	colsum = sum(A(:,j)) + inc;
	
	A(:,j) /= colsum;

	bb = ( floor( (j + b - 1) / b ) - 1 ) * b + 1; 
	be = bb + b - 1;

	A(j,bb:j-1) /= colsum;
	A(j,j+1:be) /= colsum;
end
