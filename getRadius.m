function radius = getRadius(A) 
N = size(A,1);
rv = ones(N,1)/sqrt(N);
for riter = 1:1000
   pre_rv = rv;
   rv = A * rv;
   radius = sqrt(rv'*rv);
   rv = rv./radius;
   err = min(sqrt(sum((rv-pre_rv).^2)),sqrt(sum((rv+pre_rv).^2)));
   % fprintf('error: %f, radius: %f\n', err, radius);
  if err < 0.00001
     break;
  end
end
end