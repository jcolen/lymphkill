function [] = dose_script(maskname, filesname, time, gated, outfile)

load(maskname);

[blood, blood_1frac] = LymphKill_masks(masks, time, gated, filesname);

save(outfile, 'blood', 'blood_1frac');

end
