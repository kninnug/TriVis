import fs from 'fs';
import tape from 'tape';
import {orient2d} from 'robust-predicates';
import Delaunator from 'delaunator';
import Constrainautor from '@kninnug/constrainautor';
import triangularExpansion from './TriVis.mjs';

const EPSILON = 2**-50;
const references = [{
	file: './tests/amitexp.json',
	query: [563, 241],
	output: [
		[770, 270, 590, 270],
		[780, 188.03149606299212, 780, 271.4009661835749],
		[670, 210, 690, 210],
		[530, 210, 670, 210],
		[510, 210, 530, 210],
		[470, 190, 470, 186.60377358490567],
		[450, 190, 470, 190],
		[446, 190, 450, 190],
		[434, 190, 446, 190],
		[430, 190, 434, 190],
		[410, 190, 430, 190],
		[370, 270, 370, 176.66666666666669],
		[450, 270, 370, 270],
		[450, 290, 450, 270],
		[370, 346.60377358490564, 370, 324.69026548672565],
		[530, 270, 510, 270],
		[530, 430, 530, 270],
		[560, 430, 530, 430],
		[560, 450, 560, 430],
		[600.0745610164117, 780.0000000000011, 555.2631578947369, 779.9999999999975],
		[591.3911190204142, 590, 588.233031891709, 607.8446454055273],
		[597.3508893593265, 556.3245553203368, 591.3911190204142, 590],
		[602.570792255239, 538.3161194238622, 597.3508893593265, 556.3245553203368],
		[606.5835921350013, 524.4721359549995, 602.570792255239, 538.3161194238622],
		[614.9343910395271, 507.7705381459478, 606.5835921350013, 524.4721359549995],
		[617.0585914563812, 503.52213731223964, 614.9343910395271, 507.7705381459478],
		[604.4721359549997, 346.58359213500125, 586.5835921350013, 355.52786404500046],
		[646.180339887499, 430, 604.4721359549997, 346.58359213500125],
		[650, 430, 646.180339887499, 430],
		[780, 634.8148148148148, 780, 712.4137931034484],
		[590, 270, 590, 290]
	]
}, {
	file: './tests/amitexp.json',
	query: [530, 234],
	output: [
		[770, 270, 590, 270],
		[780, 196.5, 780, 271.5],
		[670, 210, 690, 210],
		[530, 210, 670, 210],
		[510, 210, 530, 210],
		[470, 170, 470, 162],
		[470, 190, 470, 170],
		[450, 190, 470, 190],
		[446, 190, 450, 190],
		[434, 190, 446, 190],
		[430, 190, 434, 190],
		[410, 190, 430, 190],
		[370, 270, 370, 175.33333333333334],
		[450, 270, 370, 270],
		[450, 290, 450, 270],
		[370, 450, 370, 346],
		[376.6666666666667, 510, 325.55555555555554, 510],
		[530, 270, 510, 270],
		[560, 430, 530, 430],
		[588.6075312369494, 612.5403531267184, 589.3339210354185, 621.6482840980675],
		[588.233031891709, 607.8446454055273, 588.6075312369494, 612.5403531267184],
		[591.3911190204142, 590, 588.233031891709, 607.8446454055273],
		[597.3508893593265, 556.3245553203368, 591.3911190204142, 590],
		[602.570792255239, 538.3161194238622, 597.3508893593265, 556.3245553203368],
		[606.5835921350013, 524.4721359549995, 602.570792255239, 538.3161194238622],
		[614.9343910395271, 507.7705381459478, 606.5835921350013, 524.4721359549995],
		[630, 477.63932022500194, 614.9343910395271, 507.7705381459478],
		[630, 448.7758024182177, 630, 477.63932022500194],
		[604.4721359549997, 346.58359213500125, 586.5835921350013, 355.52786404500046],
		[780, 467.33333333333337, 780, 611.9386433975478],
		[590, 270, 590, 290]
	]
}, {
	file: './tests/amitexp.json',
	query: [387, 479],
	output: [
		[430, 510, 290, 510],
		[430, 470, 430, 510],
		[510, 450, 525.5555555555555, 450],
		[510, 430, 510, 450],
		[510, 270, 510, 430],
		[530, 210, 545.311004784689, 210],
		[510, 210, 530, 210],
		[510, 190, 510, 210],
		[510, 110, 510, 190],
		[370, 290, 450, 290],
		[370, 450, 370, 290],
		[350, 450, 370, 450],
		[193.97028311431427, 345.57988645363065, 179.27433670388717, 316.1879936327764],
		[199.7103162846646, 357.0599527943313, 198.6824314212446, 348.83687388697126],
		[203.41640786499872, 364.47213595499954, 199.7103162846646, 357.0599527943313],
		[200.800343220057, 365.7801682774704, 203.41640786499872, 364.47213595499954],
		[203.7592105464341, 386.43506512139453, 200.800343220057, 365.7801682774704],
		[290, 510, 290, 430]
	]
}, {
	file: './tests/amitexp.json',
	query: [530, 233, 740.0107140124046, 496.75082938258225, 266.2491706174178, 443.0107140124046],
	output: [
		[370, 450, 370, 360.3994638069705],
		[380.27027027027026, 510, 325.7603686635945, 510],
		[530, 270, 510, 270],
		[560, 430, 530, 430],
		[588.6075312369494, 612.5403531267184, 589.4976524169139, 623.7012508710677],
		[588.233031891709, 607.8446454055273, 588.6075312369494, 612.5403531267184],
		[591.3911190204142, 590, 588.233031891709, 607.8446454055273],
		[597.3508893593265, 556.3245553203368, 591.3911190204142, 590],
		[602.570792255239, 538.3161194238622, 597.3508893593265, 556.3245553203368],
		[606.5835921350013, 524.4721359549995, 602.570792255239, 538.3161194238622],
		[614.9343910395271, 507.7705381459478, 606.5835921350013, 524.4721359549995],
		[630, 477.63932022500194, 614.9343910395271, 507.7705381459478],
		[630, 449.5430991950185, 630, 477.63932022500194],
		[604.4721359549997, 346.58359213500125, 586.5835921350013, 355.52786404500046],
		[780, 546.9730639730642, 780, 614.2956036457548]
	]
}, {
	file: './tests/amitp.json',
	query: [550, 536],
	output: [
		[620, 520, 610, 550],
		[640, 480, 620, 520],
		[640, 440, 640, 480],
		[600, 360, 640, 440],
		[600, 280, 622.7272727272727, 279.9999999999999],
		[550, 200, 615.625, 200],
		[520, 440, 550, 440],
		[423.19101123595505, 140, 426.25, 140],
		[420, 180, 436, 180],
		[323.984375, 10, 357.92134831460675, 10],
		[360, 280, 440, 280],
		[360, 439.27272727272725, 360, 280],
		[440, 520, 440, 480],
		[440, 536, 440, 520],
		[280, 545.8181818181819, 280, 536],
		[440, 556, 440, 540],
		[280, 594.9090909090909, 280, 585.0909090909091],
		[440, 576, 440, 560],
		[280, 644, 280, 634.1818181818182],
		[440, 596, 440, 580],
		[113.4375, 790, 84.33333333333331, 790],
		[440, 616, 440, 600],
		[217.3809523809524, 790, 200.75, 790],
		[440, 636, 440, 620],
		[281.34615384615387, 790, 270.6, 790],
		[440, 656, 440, 640],
		[324.67741935483866, 790, 317.16666666666663, 790],
		[440, 676, 440, 660],
		[355.9722222222223, 790, 350.42857142857144, 790],
		[440, 760, 440, 680],
		[672.9032258064516, 790, 425.26785714285717, 790],
		[600, 620, 610, 660],
		[600, 600, 600, 620],
		[610, 550, 600, 600]
	]
}, {
	file: './tests/amitp.json',
	query: [550, 536, 479.28932188134524, 536, 550, 465.28932188134524],
	output: [
		[520, 440, 550, 440],
		[423.19101123595505, 140, 426.25, 140],
		[420, 180, 436, 180],
		[323.984375, 10, 357.92134831460675, 10],
		[360, 280, 440, 280],
		[360, 439.27272727272725, 360, 280],
		[440, 520, 440, 480],
		[440, 536, 440, 520]
	]
}, {
	file: './tests/amitp.json',
	query: [332, 386, 380.7903679018718, 290.5405845398161, 427.4594154601839, 434.7903679018718],
	output: [[360, 331.2173913043478, 360, 400.31111111111113]]
}, {
	file: './tests/amitp.json',
	query: [332, 386, 345.6256212499073, 264.21600086401304, 423.0689887208956, 304.0034190142456],
	output: [
		[360, 280, 360, 360.7893953820272],
		[343.8596520293065, 280.0000000000001, 360, 280]
	]
}, {
	file: './tests/ipa/mei-2.json',
	query: [459, 440],
	output: [
		[332, 462, 450, 388],
		[347, 536, 332, 462],
		[386, 639, 254.27586206896552, 615.4778325123152],
		[459, 602.8724489795918, 386, 639],
		[533, 411, 459, 514],
		[730.3329161451815, 331.46683354192743, 731.9184706332667, 333.04546421128737],
		[561, 313, 609, 380],
		[532, 134, 627.9989491583714, 229.57973977340035],
		[416.73751474636254, 195.81675186787243, 532, 134]
	]
}, {
	file: './tests/strain.json',
	query: [97, 156],
	output: [
		[5, 201, 53, 98],
		[194, 288, 5, 201],
		[294.9739572736521, 216.60427263479147, 194, 288],
		[413, 43, 146, 171],
		[278, 5, 413, 43],
		[196.30195510499638, 38.76852522326817, 278, 5],
		[53, 98, 196.30195510499638, 38.76852522326817]
	]
}, {
	file: './tests/strain.json',
	query: [330, 109],
	output: [
		[392, 148, 248.2189868368568, 249.66334264060632],
		[413, 43, 392, 148],
		[146, 171, 413, 43],
		[194, 288, 26.963427829474142, 211.11014931832938],
		[248.2189868368568, 249.66334264060632, 194, 288]
	]
}, {
	file: './tests/strain.json',
	query: [102, 192],
	output: [
		[194, 288, 5, 201],
		[324.46860547150584, 195.74947087873323, 194, 288],
		[392, 148, 324.46860547150584, 195.74947087873323],
		[412.8743718592965, 43.62814070351757, 392, 148],
		[413, 43, 146, 171],
		[278, 5, 413, 43],
		[184.0410117176336, 43.836381823378105, 278, 5],
		[53, 98, 184.0410117176336, 43.836381823378105],
		[5, 201, 53, 98]
	]
}, {
	file: './tests/strain.json',
	query: [146, 171],
	output: []
}, {
	file: './tests/tri.json',
	query: [125, 118],
	output: []
}, {
	file: './tests/ipa/matisse-nuit.json',
	query: [457, 249],
	output: [
		[453, 234, 460, 250],
		[401, 138, 453, 234],
		[371, 99, 375.6817975487971, 87.81570585565134],
		[344.61708394698087, 89.79086892488957, 371, 99],
		[445, 238, 445, 232],
		[440, 242, 445, 238],
		[421, 275, 440, 242],
		[396, 359, 421, 275],
		[377, 415, 396, 359],
		[362, 472, 377, 415],
		[345, 547, 362, 472],
		[299, 689, 293.8902272105119, 682.9885026006023],
		[309.54466230936816, 701.653594771242, 299, 689],
		[337, 637, 328, 645],
		[337.19060585432265, 637.0190605854323, 337, 637],
		[383, 492, 369, 534],
		[401, 446, 383, 492],
		[420, 398, 401, 446],
		[441, 333, 420, 398],
		[455, 281, 441, 333],
		[460, 250, 455, 281]
	]
}, {
	file: './tests/ipa/matisse-nuit.json',
	query: [352, 602],
	output: [
		[318, 610, 323, 595],
		[271, 623, 271.23023255813956, 621.0046511627907],
		[282, 669, 271, 623],
		[299, 689, 282, 669],
		[301.6601671309192, 692.1922005571031, 299, 689],
		[337, 637, 328, 645],
		[347, 638, 337, 637],
		[355, 647, 347, 638],
		[364, 685, 355, 647],
		[377, 759, 375.2477064220183, 762.7966360856267],
		[389, 749, 377, 759],
		[399, 762, 389, 749],
		[418, 768, 399, 762],
		[430.3447136563877, 757.3386563876652, 418, 768],
		[405, 704, 410, 717],
		[407, 684, 405, 704],
		[412.9230769230769, 672.1538461538462, 407, 684],
		[379, 633, 385, 640],
		[380, 614, 379, 633],
		[395, 576, 380, 614],
		[428, 509, 395, 576],
		[455, 461, 428, 509],
		[500, 378, 455, 461],
		[510.28985507246375, 355.36231884057975, 500, 378],
		[394.2061855670103, 535.9381443298969, 395, 535],
		[365, 551, 375, 566],
		[369, 534, 365, 551],
		[440, 242, 442.5, 240],
		[421, 275, 440, 242],
		[396, 359, 421, 275],
		[377, 415, 396, 359],
		[362, 472, 377, 415],
		[345, 547, 362, 472],
		[336, 550, 345, 547],
		[323, 527, 336, 550],
		[302, 520, 305.3826429980276, 481.43786982248525],
		[319, 571, 302, 520],
		[323, 595, 319, 571]
	]
}, {
	file: './tests/ipa/seidel-3.json',
	query: [256, 352],
	output: [
		[175, 282, 247, 315],
		[112, 746, 114.77940300483155, 229.9575087696075],
		[198, 657, 112, 746],
		[285, 780, 198, 657],
		[294, 490, 285, 780],
		[384, 583, 294, 490],
		[430, 650, 384, 583],
		[513, 540, 430, 650],
		[616, 516, 582.9042881165919, 591.1362107623319],
		[707, 287, 616, 516],
		[579, 216, 707, 287],
		[594, 101, 579, 216],
		[506, 155, 594, 101],
		[505, 17, 506, 155],
		[451.89786812097174, 70.79176995537927, 505, 17],
		[247, 315, 318, 263]
	]
}, {
	file: './tests/ipa/seidel-3.json',
	query: [256, 352, 413.2187408195795, 485.0122833993383, 234.99260397048673, 556.8626108201768],
	output: [
		[285, 780, 221.34034829332205, 689.9984234491795],
		[294, 490, 285, 780],
		[384, 583, 294, 490],
		[430, 650, 384, 583],
		[499.44587113039216, 557.9633033211669, 430, 650]
	]
}, {
	file: './tests/ipa/seidel-3.json',
	query: [256, 352, 355.5263016439092, 94.86284733417926, 524.4212271396464, 288.9520434839767],
	output: [
		[579, 216, 655.1581590242263, 258.2439788337506],
		[594, 101, 579, 216],
		[506, 155, 594, 101],
		[505, 17, 506, 155],
		[451.89786812097174, 70.79176995537927, 505, 17],
		[279.5475259828696, 291.16237533648984, 318, 263]
	]
}, {
	file: './tests/ipa/seidel-3.json',
	query: [256, 352, 182.2758425261011, 355.8403912253514, 211.2025977215036, 293.3212751578686],
	output: [
		[175, 282, 217.393275160481, 301.43025111522047],
		[114.08226559954085, 359.3926870185832, 114.77940300483158, 229.95750876960753]
	]
}];

function sqdist(x1, y1, x2, y2){
	const dx = x2 - x1,
		dy = y2 - y1;
	
	return dx * dx + dy * dy;
}

function validatePoly(t, qx, qy, poly){
	const len = poly.length;
	if(!len){
		return;
	}
	
	let [plx, ply, prx, pry] = poly[0];
	for(let i = 1; i < len; i++){
		const [lx, ly, rx, ry] = poly[i],
			oSeg = orient2d(qx, qy, rx, ry, lx, ly),
			oPrev = orient2d(qx, qy, prx, pry, rx, ry);
		if(oSeg < 0){
			t.fail(`segment ${i} [${lx}, ${ly}, ${rx}, ${ry}] not oriented left-right (${oSeg})`);
		}
		if(oPrev < 0){
			t.fail(`segment ${i}: [${lx}, ${ly}, ${rx}, ${ry}] not left of previous: [${plx}, ${ply}, ${prx}, ${pry}] (${oPrev})`);
		}
		([plx, ply, prx, pry] = [lx, ly, rx, ry]);
	}
	t.pass(`segments oriented counter-clockwise`);
}

function comparePolys(t, poly, ref){
	if(poly.length !== ref.length){
		t.fail(`wrong number of segments: ${poly.length} !== ${ref.length}`);
		console.log('poly', poly);
		return;
	}
	
	const len = poly.length;
	let maxErr = 0;
	for(let i = 0; i < len; i++){
		const p = poly[i],
			r = ref[i],
			dl = sqdist(p[0], p[1], r[0], r[1]),
			dr = sqdist(p[2], p[3], r[2], r[3]);
		if(dl > EPSILON || dr > EPSILON){
			t.fail(`segment ${i} not equal to reference`);
		}
		maxErr = Math.max(maxErr, dl, dr);
	}
	t.pass(`visiblity polygon equal to reference, maximum error: ${maxErr}`);
}

function testFile(t, json, name){
	const {points, edges} = json,
		del = Delaunator.from(points),
		con = new Constrainautor(del);
	con.constrainAll(edges);
	
	for(const {file, query: [qx, qy, ilx = NaN, ily = NaN, irx = NaN, iry = NaN], output} of references){
		if(file !== name){
			continue;
		}
		
		const vis = triangularExpansion(del, qx, qy, 
				edg => con.isConstrained(edg), ilx, ily, irx, iry);
		
		validatePoly(t, qx, qy, vis);
		comparePolys(t, vis, output);
	}
	
	t.end();
}

const files = fs.readdirSync('./tests/', 'utf8').map(f => './tests/' + f)
		.concat(fs.readdirSync('./tests/ipa/', 'utf8').map(f => './tests/ipa/' + f))
		.filter(f => f.endsWith('.json'));

function main(args){
	//if(!args.length){
	//	tape.test("Example", testExample);
	//}
	
	args = args.length ? args : files;
	
	for(const file of args){
		const json = JSON.parse(fs.readFileSync(file, 'utf8'));
		tape.test(file, (t) => testFile(t, json, file));
	}
}

main(process.argv.slice(2));
