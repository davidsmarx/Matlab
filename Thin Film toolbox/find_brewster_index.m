function n = find_brewster_index(ns,theta)

n = fzero(@reflection,ns,[],ns,theta);

return


function r = reflection(n,ns,theta)

[R, T, r, tt] = thin_film_filter_2([ns n],[],theta,1,1); % independent of lam since d = []

return