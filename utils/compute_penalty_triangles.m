function penalty_triang = compute_penalty_triangles(B1,B02)
%Created by AGM
%Function to check the cost (in terms of links) of activating a triangle. If the 3 links are present, the cost is zero
%
%Output: Vector with lengh maximum number of triangles (N choose 3)
%Inputs: Full incidence matrices (B1: NxE current estimate)  (B02: ExT all possible triangles) 

plot_boolean = 0;

NT=size(B02,2);
w_edge = sum(abs(B1))./sum(sum(abs(B1)));
w_triang=zeros(1,NT);
for t=1:NT
    w_triang(t) = sum(abs(B02(:,t).*w_edge'));    
end



w_triang2=w_triang+min(w_triang(w_triang>0))/3;
penalty_triang=1./w_triang2-1/max(w_triang2);

if plot_boolean
    figure
    subplot(211)
    stem(w_triang)
    ylabel('Clique likelihood')
    subplot(212)
    stem(penalty_triang)
    ylabel('Cost creating clique')
end

end