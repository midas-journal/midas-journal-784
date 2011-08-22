REM Test using down-sized CT
getDRRSiddonJacobsRayTracing -rp 0 -rx -3 -ry 4 -rz 2 -t 5 5 5 -iso 99.62 101.18 65 -res 1 1 -size 256 256 -o boxheadDRRDev1_G0.tif BoxheadCT.img
getDRRSiddonJacobsRayTracing -rp 90 -rx -3 -ry 4 -rz 2 -t 5 5 5 -iso 99.62 101.18 65 -res 1 1 -size 256 256 -o boxheadDRRDev1_G90.tif BoxheadCT.img
2D3DTwoProjRegistration -res 1 1 1 1 -iso 99.62 101.18 65 -o boxheadDRRDev1_G0_Reg.tif boxheadDRRDev1_G90_Reg.tif boxheadDRRDev1_G0.tif 0 boxheadDRRDev1_G90.tif 90 BoxheadCT.img

REM Test using full size CT
getDRRSiddonJacobsRayTracing -rp 0 -rx -3 -ry 4 -rz 2 -t 5 5 5 -iso 255 259 130 -res 0.5 0.5 -size 512 512 -o BoxheadDRRFullDev1_G0.tif BoxheadCTFull.img
getDRRSiddonJacobsRayTracing -rp 90 -rx -3 -ry 4 -rz 2 -t 5 5 5 -iso 255 259 130 -res 0.5 0.5 -size 512 512 -o BoxheadDRRFullDev1_G90.tif BoxheadCTFull.img
2D3DTwoProjRegistration -res 0.5 0.5 0.5 0.5 -iso 255 259 130 -o BoxheadDRRFullDev1_G0_Reg.tif BoxheadDRRFullDev1_G90_Reg.tif BoxheadDRRFullDev1_G0.tif 0 BoxheadDRRFullDev1_G90.tif 90 BoxheadCTFull.img
