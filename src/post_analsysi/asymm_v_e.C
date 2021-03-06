{
//=========Macro generated from canvas: cGreat/Energy vs. A
//=========  (Mon Aug  3 13:42:42 2009) by ROOT version5.18/00
   TCanvas *cGreat = new TCanvas("cGreat", "Energy vs. A",510,97,700,502);
   gStyle->SetOptStat(0);
   cGreat->SetHighLightColor(2);
   cGreat->Range(-130.347,-0.07847063,654.2318,0.1573557);
   cGreat->SetFillColor(0);
   cGreat->SetBorderMode(0);
   cGreat->SetBorderSize(2);
   cGreat->SetGridy();
   cGreat->SetLeftMargin(0.1393678);
   cGreat->SetFrameBorderMode(0);
   cGreat->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(160);
   gre->SetName("Graph");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineColor(2);
   gre->SetMarkerColor(2);
   gre->SetMarkerStyle(20);
   gre->SetPoint(0,-5,0.006389111);
   gre->SetPointError(0,0.1,0.04563862);
   gre->SetPoint(1,5,nan);
   gre->SetPointError(1,0.1,nan);
   gre->SetPoint(2,13.94332,0.06942048);
   gre->SetPointError(2,0.1,0.04703041);
   gre->SetPoint(3,24.08874,0.004468302);
   gre->SetPointError(3,0.1,0.01874078);
   gre->SetPoint(4,35,0.04689694);
   gre->SetPointError(4,0.1,0.01391268);
   gre->SetPoint(5,45,0.02290334);
   gre->SetPointError(5,0.1,0.01207032);
   gre->SetPoint(6,55,0.03184075);
   gre->SetPointError(6,0.1,0.01116758);
   gre->SetPoint(7,65,0.02544128);
   gre->SetPointError(7,0.1,0.01048324);
   gre->SetPoint(8,75,0.04229715);
   gre->SetPointError(8,0.1,0.01014151);
   gre->SetPoint(9,85,0.04486172);
   gre->SetPointError(9,0.1,0.009795184);
   gre->SetPoint(10,95,0.0365765);
   gre->SetPointError(10,0.1,0.009548549);
   gre->SetPoint(11,105,0.02632267);
   gre->SetPointError(11,0.1,0.009319111);
   gre->SetPoint(12,115,0.04049164);
   gre->SetPointError(12,0.1,0.009150744);
   gre->SetPoint(13,125,0.02253341);
   gre->SetPointError(13,0.1,0.009042712);
   gre->SetPoint(14,135,0.03496574);
   gre->SetPointError(14,0.1,0.008879167);
   gre->SetPoint(15,145,0.04187489);
   gre->SetPointError(15,0.1,0.008763972);
   gre->SetPoint(16,155,0.03605345);
   gre->SetPointError(16,0.1,0.008655715);
   gre->SetPoint(17,165,0.03708767);
   gre->SetPointError(17,0.1,0.008550943);
   gre->SetPoint(18,175,0.0386679);
   gre->SetPointError(18,0.1,0.008454361);
   gre->SetPoint(19,185,0.04572165);
   gre->SetPointError(19,0.1,0.008445212);
   gre->SetPoint(20,195,0.0295763);
   gre->SetPointError(20,0.1,0.008290489);
   gre->SetPoint(21,205,0.03604979);
   gre->SetPointError(21,0.1,0.008332276);
   gre->SetPoint(22,215,0.03972457);
   gre->SetPointError(22,0.1,0.00822089);
   gre->SetPoint(23,225,0.04828961);
   gre->SetPointError(23,0.1,0.008303013);
   gre->SetPoint(24,235,0.04088816);
   gre->SetPointError(24,0.1,0.008228064);
   gre->SetPoint(25,245,0.03304505);
   gre->SetPointError(25,0.1,0.008203362);
   gre->SetPoint(26,255,0.04529746);
   gre->SetPointError(26,0.1,0.008081101);
   gre->SetPoint(27,265,0.04270395);
   gre->SetPointError(27,0.1,0.008187038);
   gre->SetPoint(28,275,0.04155198);
   gre->SetPointError(28,0.1,0.008267311);
   gre->SetPoint(29,285,0.0595088);
   gre->SetPointError(29,0.1,0.008253983);
   gre->SetPoint(30,295,0.03813331);
   gre->SetPointError(30,0.1,0.008235696);
   gre->SetPoint(31,305,0.0456729);
   gre->SetPointError(31,0.1,0.008329624);
   gre->SetPoint(32,315,0.04409898);
   gre->SetPointError(32,0.1,0.008381896);
   gre->SetPoint(33,325,0.03963665);
   gre->SetPointError(33,0.1,0.008449449);
   gre->SetPoint(34,335,0.04614977);
   gre->SetPointError(34,0.1,0.008492964);
   gre->SetPoint(35,345,0.03580619);
   gre->SetPointError(35,0.1,0.00857545);
   gre->SetPoint(36,355,0.05289987);
   gre->SetPointError(36,0.1,0.008611034);
   gre->SetPoint(37,365,0.03645456);
   gre->SetPointError(37,0.1,0.008788171);
   gre->SetPoint(38,375,0.03326673);
   gre->SetPointError(38,0.1,0.008991885);
   gre->SetPoint(39,385,0.04136398);
   gre->SetPointError(39,0.1,0.009035574);
   gre->SetPoint(40,395,0.0478886);
   gre->SetPointError(40,0.1,0.009272826);
   gre->SetPoint(41,405,0.06589146);
   gre->SetPointError(41,0.1,0.009432013);
   gre->SetPoint(42,415,0.04179408);
   gre->SetPointError(42,0.1,0.009677216);
   gre->SetPoint(43,425,0.04831837);
   gre->SetPointError(43,0.1,0.009903408);
   gre->SetPoint(44,435,0.05964095);
   gre->SetPointError(44,0.1,0.01017682);
   gre->SetPoint(45,445,0.05398298);
   gre->SetPointError(45,0.1,0.0103556);
   gre->SetPoint(46,455,0.046765);
   gre->SetPointError(46,0.1,0.0105618);
   gre->SetPoint(47,465,0.03616569);
   gre->SetPointError(47,0.1,0.01120012);
   gre->SetPoint(48,475,0.05688507);
   gre->SetPointError(48,0.1,0.01139146);
   gre->SetPoint(49,485,0.05697509);
   gre->SetPointError(49,0.1,0.0118835);
   gre->SetPoint(50,495,0.05077604);
   gre->SetPointError(50,0.1,0.01246748);
   gre->SetPoint(51,505,0.03578008);
   gre->SetPointError(51,0.1,0.01285779);
   gre->SetPoint(52,515,0.06316115);
   gre->SetPointError(52,0.1,0.01340473);
   gre->SetPoint(53,525,0.05336101);
   gre->SetPointError(53,0.1,0.01429867);
   gre->SetPoint(54,535,0.03801421);
   gre->SetPointError(54,0.1,0.01492184);
   gre->SetPoint(55,545,0.04481478);
   gre->SetPointError(55,0.1,0.01629204);
   gre->SetPoint(56,555,0.09902191);
   gre->SetPointError(56,0.1,0.01736529);
   gre->SetPoint(57,565,0.07795619);
   gre->SetPointError(57,0.1,0.01873053);
   gre->SetPoint(58,575,0.06518294);
   gre->SetPointError(58,0.1,0.01967685);
   gre->SetPoint(59,585,0.07361266);
   gre->SetPointError(59,0.1,0.02192087);
   gre->SetPoint(60,595,0.094161);
   gre->SetPointError(60,0.1,0.02397361);
   gre->SetPoint(61,605,0.08118724);
   gre->SetPointError(61,0.1,0.02739183);
   gre->SetPoint(62,615,0.05681355);
   gre->SetPointError(62,0.1,0.03050748);
   gre->SetPoint(63,625,nan);
   gre->SetPointError(63,0.1,nan);
   gre->SetPoint(64,635,nan);
   gre->SetPointError(64,0.1,nan);
   gre->SetPoint(65,645,nan);
   gre->SetPointError(65,0.1,nan);
   gre->SetPoint(66,655,nan);
   gre->SetPointError(66,0.1,nan);
   gre->SetPoint(67,665,nan);
   gre->SetPointError(67,0.1,nan);
   gre->SetPoint(68,675,nan);
   gre->SetPointError(68,0.1,nan);
   gre->SetPoint(69,685,nan);
   gre->SetPointError(69,0.1,nan);
   gre->SetPoint(70,695,nan);
   gre->SetPointError(70,0.1,nan);
   gre->SetPoint(71,705,nan);
   gre->SetPointError(71,0.1,nan);
   gre->SetPoint(72,715,nan);
   gre->SetPointError(72,0.1,nan);
   gre->SetPoint(73,725,nan);
   gre->SetPointError(73,0.1,nan);
   gre->SetPoint(74,735,nan);
   gre->SetPointError(74,0.1,nan);
   gre->SetPoint(75,745,nan);
   gre->SetPointError(75,0.1,nan);
   gre->SetPoint(76,755,nan);
   gre->SetPointError(76,0.1,nan);
   gre->SetPoint(77,765,nan);
   gre->SetPointError(77,0.1,nan);
   gre->SetPoint(78,775,nan);
   gre->SetPointError(78,0.1,nan);
   gre->SetPoint(79,785,nan);
   gre->SetPointError(79,0.1,nan);
   gre->SetPoint(80,795,nan);
   gre->SetPointError(80,0.1,nan);
   gre->SetPoint(81,805,nan);
   gre->SetPointError(81,0.1,nan);
   gre->SetPoint(82,815,nan);
   gre->SetPointError(82,0.1,nan);
   gre->SetPoint(83,825,nan);
   gre->SetPointError(83,0.1,nan);
   gre->SetPoint(84,835,nan);
   gre->SetPointError(84,0.1,nan);
   gre->SetPoint(85,845,nan);
   gre->SetPointError(85,0.1,nan);
   gre->SetPoint(86,855,nan);
   gre->SetPointError(86,0.1,nan);
   gre->SetPoint(87,865,nan);
   gre->SetPointError(87,0.1,nan);
   gre->SetPoint(88,875,nan);
   gre->SetPointError(88,0.1,nan);
   gre->SetPoint(89,885,nan);
   gre->SetPointError(89,0.1,nan);
   gre->SetPoint(90,895,nan);
   gre->SetPointError(90,0.1,nan);
   gre->SetPoint(91,905,nan);
   gre->SetPointError(91,0.1,nan);
   gre->SetPoint(92,915,nan);
   gre->SetPointError(92,0.1,nan);
   gre->SetPoint(93,925,nan);
   gre->SetPointError(93,0.1,nan);
   gre->SetPoint(94,935,nan);
   gre->SetPointError(94,0.1,nan);
   gre->SetPoint(95,945,nan);
   gre->SetPointError(95,0.1,nan);
   gre->SetPoint(96,955,nan);
   gre->SetPointError(96,0.1,nan);
   gre->SetPoint(97,965,nan);
   gre->SetPointError(97,0.1,nan);
   gre->SetPoint(98,975,nan);
   gre->SetPointError(98,0.1,nan);
   gre->SetPoint(99,985,nan);
   gre->SetPointError(99,0.1,nan);
   gre->SetPoint(100,995,nan);
   gre->SetPointError(100,0.1,nan);
   gre->SetPoint(101,1005,nan);
   gre->SetPointError(101,0.1,nan);
   gre->SetPoint(102,1015,nan);
   gre->SetPointError(102,0.1,nan);
   gre->SetPoint(103,1025,nan);
   gre->SetPointError(103,0.1,nan);
   gre->SetPoint(104,1035,nan);
   gre->SetPointError(104,0.1,nan);
   gre->SetPoint(105,1045,nan);
   gre->SetPointError(105,0.1,nan);
   gre->SetPoint(106,1055,nan);
   gre->SetPointError(106,0.1,nan);
   gre->SetPoint(107,1065,nan);
   gre->SetPointError(107,0.1,nan);
   gre->SetPoint(108,1075,nan);
   gre->SetPointError(108,0.1,nan);
   gre->SetPoint(109,1085,nan);
   gre->SetPointError(109,0.1,nan);
   gre->SetPoint(110,1095,nan);
   gre->SetPointError(110,0.1,nan);
   gre->SetPoint(111,1105,nan);
   gre->SetPointError(111,0.1,nan);
   gre->SetPoint(112,1115,nan);
   gre->SetPointError(112,0.1,nan);
   gre->SetPoint(113,1125,nan);
   gre->SetPointError(113,0.1,nan);
   gre->SetPoint(114,1135,nan);
   gre->SetPointError(114,0.1,nan);
   gre->SetPoint(115,1145,nan);
   gre->SetPointError(115,0.1,nan);
   gre->SetPoint(116,1155,nan);
   gre->SetPointError(116,0.1,nan);
   gre->SetPoint(117,1165,nan);
   gre->SetPointError(117,0.1,nan);
   gre->SetPoint(118,1175,nan);
   gre->SetPointError(118,0.1,nan);
   gre->SetPoint(119,1185,nan);
   gre->SetPointError(119,0.1,nan);
   gre->SetPoint(120,1195,nan);
   gre->SetPointError(120,0.1,nan);
   gre->SetPoint(121,1205,nan);
   gre->SetPointError(121,0.1,nan);
   gre->SetPoint(122,1215,nan);
   gre->SetPointError(122,0.1,nan);
   gre->SetPoint(123,1225,nan);
   gre->SetPointError(123,0.1,nan);
   gre->SetPoint(124,1235,nan);
   gre->SetPointError(124,0.1,nan);
   gre->SetPoint(125,1245,nan);
   gre->SetPointError(125,0.1,nan);
   gre->SetPoint(126,1255,nan);
   gre->SetPointError(126,0.1,nan);
   gre->SetPoint(127,1265,nan);
   gre->SetPointError(127,0.1,nan);
   gre->SetPoint(128,1275,nan);
   gre->SetPointError(128,0.1,nan);
   gre->SetPoint(129,1285,nan);
   gre->SetPointError(129,0.1,nan);
   gre->SetPoint(130,1295,nan);
   gre->SetPointError(130,0.1,nan);
   gre->SetPoint(131,1305,nan);
   gre->SetPointError(131,0.1,nan);
   gre->SetPoint(132,1315,nan);
   gre->SetPointError(132,0.1,nan);
   gre->SetPoint(133,1325,nan);
   gre->SetPointError(133,0.1,nan);
   gre->SetPoint(134,1335,nan);
   gre->SetPointError(134,0.1,nan);
   gre->SetPoint(135,1345,nan);
   gre->SetPointError(135,0.1,nan);
   gre->SetPoint(136,1355,nan);
   gre->SetPointError(136,0.1,nan);
   gre->SetPoint(137,1365,nan);
   gre->SetPointError(137,0.1,nan);
   gre->SetPoint(138,1375,nan);
   gre->SetPointError(138,0.1,nan);
   gre->SetPoint(139,1385,nan);
   gre->SetPointError(139,0.1,nan);
   gre->SetPoint(140,1395,nan);
   gre->SetPointError(140,0.1,nan);
   gre->SetPoint(141,1405,nan);
   gre->SetPointError(141,0.1,nan);
   gre->SetPoint(142,1415,nan);
   gre->SetPointError(142,0.1,nan);
   gre->SetPoint(143,1425,nan);
   gre->SetPointError(143,0.1,nan);
   gre->SetPoint(144,1435,nan);
   gre->SetPointError(144,0.1,nan);
   gre->SetPoint(145,1445,nan);
   gre->SetPointError(145,0.1,nan);
   gre->SetPoint(146,1455,nan);
   gre->SetPointError(146,0.1,nan);
   gre->SetPoint(147,1465,nan);
   gre->SetPointError(147,0.1,nan);
   gre->SetPoint(148,1475,nan);
   gre->SetPointError(148,0.1,nan);
   gre->SetPoint(149,1485,nan);
   gre->SetPointError(149,0.1,nan);
   gre->SetPoint(150,1495,nan);
   gre->SetPointError(150,0.1,nan);
   gre->SetPoint(151,1505,nan);
   gre->SetPointError(151,0.1,nan);
   gre->SetPoint(152,1515,nan);
   gre->SetPointError(152,0.1,nan);
   gre->SetPoint(153,1525,nan);
   gre->SetPointError(153,0.1,nan);
   gre->SetPoint(154,1535,nan);
   gre->SetPointError(154,0.1,nan);
   gre->SetPoint(155,1545,nan);
   gre->SetPointError(155,0.1,nan);
   gre->SetPoint(156,1555,nan);
   gre->SetPointError(156,0.1,nan);
   gre->SetPoint(157,1565,nan);
   gre->SetPointError(157,0.1,nan);
   gre->SetPoint(158,1575,nan);
   gre->SetPointError(158,0.1,nan);
   gre->SetPoint(159,1585,nan);
   gre->SetPointError(159,0.1,0);
   
   TH1 *Graph_h1 = new TH1F("Graph_h1","",160,-164.12,1744.12);
   Graph_h1->SetMinimum(-0.05498792);
   Graph_h1->SetMaximum(0.133873);
   Graph_h1->SetDirectory(0);
   Graph_h1->SetStats(0);
   Graph_h1->GetXaxis()->SetTitle("ADC Channel");
   Graph_h1->GetXaxis()->SetRange(13,62);
   Graph_h1->GetYaxis()->SetTitle("Raw Asymmetry ");
   Graph_h1->GetYaxis()->CenterTitle(true);
   Graph_h1->GetYaxis()->SetTitleOffset(1.6);
   gre->SetHistogram(Graph_h1);
   
   
   TF1 *fAbeta = new TF1("fAbeta","0.5*[0]*sqrt(x*(x+1022.))/(x+511.)",100,600);
   fAbeta->SetFillColor(19);
   fAbeta->SetFillStyle(0);
   fAbeta->SetLineWidth(3);
   fAbeta->SetChisquare(44.0361);
   fAbeta->GetYaxis()->SetTitleOffset(1.6);
   fAbeta->SetParameter(0,0.1154211);
   fAbeta->SetParError(0,0.003584142);
   fAbeta->SetParLimits(0,0,0);
   gre->GetListOfFunctions()->Add(fAbeta);
   gre->Draw("ap");
   cGreat->Modified();
   cGreat->cd();
   cGreat->SetSelected(cGreat);
}
