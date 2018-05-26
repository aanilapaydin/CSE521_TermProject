function showTheat(x,t,T,pflag)
% showTheat Plot and print solutions to the heat equation
%
% Synopsis: showTheat(x,t,T)
% showTheat(x,t,T,pflag)
%
% Input: x = vector of positions at which temperature is known
% t = vector of times at which solution is to be evaluated
% T = matrix of T(x,t) values. T(:,t(j)) is T(x) at time t(j)
% pflag = flag to control printing of results.
% Default: pflag = 0, do not print results
% --- Define string matrix of line styles that can be reused.
% nsymb is the total number of line style strings (plot symbols)
% In plot loop the statement isymb = 1 + rem(j-1,nsymb) is an
% index in the range 1 <= isymb <= nsymb
lineSymb = strvcat(’b-.’,’k-o’,’m-v’,’b-s’,’k*-’,’m-d’,’b-+’,’k-<’,’m-h’);
nsymb = size(lineSymb,1);
% --- Plot T(x,t): each time is a different curve with different symbol
for j=1:length(t)
	hold on
	isymb = 1 + rem(j-1,nsymb); % cyclic index for line styles.
	plot(x,T(:,j),lineSymb(isymb,:));
	% --- Build string matrix for legend. Each row is a legend string.
	s = sprintf(’t = %-4.2f’,t(j));
	if j==1
		legstr = s;
	else
		legstr = strvcat(legstr,s);
	end
end
hold off

legend(legstr,2); xlabel(’x’); ylabel(’T(x,t)’);