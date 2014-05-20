function varna(sequence, structure, filename)
VARNA = 'java -cp /home/tsuname/varna/bin/VARNA.jar fr.orsay.lri.varna.applications.VARNAcmd ';
system( [VARNA, ' -sequenceDBN ', sequence, ' -structureDBN "', structure,  '" -o ', filename]);
end