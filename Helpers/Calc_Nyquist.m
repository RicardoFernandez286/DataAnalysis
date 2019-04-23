function nyq_speed = Calc_Nyquist(reprate,WL)

nyq_speed = [num2str((WL./4)./(10^6./reprate)/2) ' mm/s (roundtrip)'];