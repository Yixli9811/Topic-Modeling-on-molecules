
import OCL from "openchemlib"; 
import fs from 'fs'; 

function getMoleculeDataFromSmiles(smiles) {

	var molecule = OCL.Molecule.fromSmiles(smiles);
	var mf = molecule.getMolecularFormula();

	var properties = new OCL.MoleculeProperties(molecule);
	var result = {
		smiles: smiles, 
		mf: mf.formula,
		mw: mf.relativeWeight,
		em: mf.absoluteWeight,
		logP: properties.logP,
		logS: properties.logS,
		psa: properties.polarSurfaceArea,
		donorCount: properties.donorCount, 
		centerCount: properties.stereoCenterCount, 
		rotatableBondCount: properties.rotatableBondCount,
		acceptorCount: properties.acceptorCount,

		};

	// return result 
	return result;
}

    
// writeFile function is defined. 

// Data which will write in a file. 

	var lines = fs.readFileSync('./cleaned_meltingPt.txt').toString('utf-8').split('\n');
	var arrayLength = lines.length;

	for (let i = 0; i < arrayLength; i++) {

    try {
	  getMoleculeDataFromSmiles(lines[i]);
	}
	catch(err) {
	  console.log(i);
	}
	}; 

	var result = lines.map(getMoleculeDataFromSmiles);

	fs.writeFile('./meltPt_prop_' + i + '.csv', convertToCSV(result), function(err) {
		if (err) {
			return console.error(err);
		}
	});

function convertToCSV(arr) {
  const array = [Object.keys(arr[0])].concat(arr)

  return array.map(it => {
    return Object.values(it).toString()
  }).join('\n')
}
