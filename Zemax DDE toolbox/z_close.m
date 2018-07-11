function z_close(zchan)
% z_close(zchan)

rc = ddeterm(zchan);

if rc == 0,
    error('Failure Terminating DDE');
end

