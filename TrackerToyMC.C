struct PADMEGeometry
{
  //coordinates are wrt the centre of the PADME dipole
  double target_z       = -732; //mm, Run 4 position (run 3 = -1028 mm, run 4 = 510 mm)
  double magnet_front_z = -490; //mm
  double magnet_back_z  = 510;  //mm
  double chamber_back_z = 2300; //mm
  double ECal_front_z   = 2436; //mm
  double ECal_length    = 230;  //mm
  double magnet_width_x = 435;  //mm, == chamber width
  double magnet_width_y = 205;  //mm, == chamber width
  double ECal_RMax      = 270;  //mm, maximum radius of highest energy cluster
  double ECal_cell_size = 21;   //mm
  int    nECal_cells    = 29;   // per row or column
  double ECal_edge      = (nECal_cells/2+0.5)*ECal_cell_size; // 304.5mm
  double deltaZ = ECal_front_z-target_z; //target-ECal distance
  double beamspotsize_X = 1.;  //mm, sigma beam X position at target
  double beamspotsize_Y = 0.5; //mm, sigma beam Y position at target
  double beamdivergence = 2;   //mrad
};

PADMEGeometry Geometry;
const Double_t me = 0.511; //MeV
const Double_t SigmaPOverP = 0.25e-2; //Run 3 paper: relative energy spread = 0.25% DOI 10.1007/JHEP11(2025)007
const Double_t fConstantMagneticField = 0.05;//0.36;//0.4542; // [T]
const Double_t position_res = 100e-6; //assume 100 microns for straws
const Int_t    nTrackerPlanes = 4;//3;
const Double_t ThetaBremPos = 14; //mrad, >98% of bremsstrahlung positrons have theta less than this at 430 MeV beam energy

const Double_t z1 = Geometry.magnet_front_z; //first scintillator position at entrance to magnet
const Double_t z2 = Geometry.chamber_back_z; //second scintillator position at back of chamber
std::vector<std::pair<TString, Double_t>> vpZPositions = {{"Z1",z1}, {"Z2",z2}};

std::vector<TString>     aTags = {"Ele","Pos","Brem", "Beam"};
std::array<TString, 2>   aXY   = {"X","Y"};

std::vector<std::pair<TString, TString>> vpsCuts;
vpsCuts.push_back({"BeamDivCut", "with hole for beam divergence"});
vpsCuts.push_back({"BremPosCut", "with hole for bremsstrahlung positron"});

const bool BeamEnergySpread = 1;
const bool ElectronEnergySpread = 1;

TRandom3 rand3 = TRandom3(1);

void FillPositionHistos(std::map<TString, TH1F*>& mh1, std::map<TString, TH2F*>& mh2, TString tag, double x, double y, TString sZKey, TString sCutName)
{
  if(sCutName.Length() == 0)
    {
      mh2.at("h" + tag + "_XY" + sZKey)->Fill(x,y);
      mh1.at("h" + tag + "_X"  + sZKey)->Fill(x);
      mh1.at("h" + tag + "_Y"  + sZKey)->Fill(y);
    }
  else 
    {
      mh2.at("h" + tag + "_XY" + sZKey + "_" + sCutName)->Fill(x,y);
      mh1.at("h" + tag + "_X"  + sZKey + "_" + sCutName)->Fill(x);
      mh1.at("h" + tag + "_Y"  + sZKey + "_" + sCutName)->Fill(y);
    }
}

//returns theta angle and energy of Bremsstrahlung positron randomly sampled from 2D histogram
//as of 1/7/26, 2D histogram used comes from MC truth for PADME run 660 conditions => ~430 MeV beam energy
std::pair<double, double> BremThetaEneSampler(TString sBremBinFile, std::map<TString, TH2F*>& mh2)
{
  int nx, ny;
  std::vector<double> xEdges, yEdges, mycdf;

  std::ifstream in(sBremBinFile, std::ios::binary);
  if (!in)
    {
      std::cerr << "Failed to open BremSampler.bin\n";
      return {0,0};
    }
  
  in.read(reinterpret_cast<char*>(&nx), sizeof(int));
  in.read(reinterpret_cast<char*>(&ny), sizeof(int));
  
  xEdges.resize(nx + 1);
  yEdges.resize(ny + 1);
  mycdf.resize(nx * ny);
  
  in.read(reinterpret_cast<char*>(xEdges.data()), xEdges.size() * sizeof(double));
  in.read(reinterpret_cast<char*>(yEdges.data()), yEdges.size() * sizeof(double));
  in.read(reinterpret_cast<char*>(mycdf.data()), mycdf.size() * sizeof(double));

  if (mycdf.empty())
    {
      std::cerr << "MYCDF is empty\n";
      return {0,0};
    }
  
  double u = rand3.Uniform(mycdf.back());
  
  int k = std::lower_bound(mycdf.begin(), mycdf.end(), u) - mycdf.begin();
  
  int ix = k % nx;
  int iy = k / nx;
  
  double theta  = rand3.Uniform(xEdges[ix], xEdges[ix+1]);
  double energy = rand3.Uniform(yEdges[iy], yEdges[iy+1]);

  mh2.at("hThetaVsEnergy")->Fill(theta,energy);
  return {theta, energy};
}

static Double_t SignalShape(double signalpeak,double bes,double lorewidth,double massn, double sqrts){ // return number of signal events produced(not the accepted) / POT at sqrts(s) at gven = 1
  const double Wb = -7E-6;// MeV
  
  double eRes = (massn*massn - 2*me*me)/(2.*(me+Wb));
  double eBeam = (sqrts*sqrts - 2*me*me)/(2.*me); // s = 2me^2 + 2meEbeam -> eBeam = (s-2me^2) / 2me
  return signalpeak*TMath::Voigt(eBeam-eRes,bes*eRes,lorewidth*2,5); // last input was 4 but with low BES can have rounding problems
}

//minimum radius of cluster for ECal acceptance given sqrts
double kineRMin(double sqrts){
  double deltaZ = Geometry.deltaZ;//target-ECal distanec
  
  double fXTarg = 0;
  double fYTarg = 0;
  //  double fZTarg = z0;//-1028;
  double fZTarg = Geometry.target_z;
  double fXCOG = 0;
  double fYCOG = 0;
  //  double fZECal = 2612.4;
  const double fMe = 0.511;
  //  double fSqrts = x[0];
  double fSqrts = sqrts;
  
  double fBeamMomentum = (fSqrts*fSqrts-2.*fMe*fMe)/(2.*fMe);
  //  double fRadiusMax = par[0]; //   fRadiusMax = 270.0;
  double fRadiusMax = Geometry.ECal_RMax;
  
  TVector3 fRTarg; fRTarg.SetXYZ(fXTarg,fYTarg,fZTarg);
  //  TVector3 fCOGAtECal; fCOGAtECal.SetXYZ(fXCOG,fYCOG,ecalZ);
  TVector3 fCOGAtECal; fCOGAtECal.SetXYZ(fXCOG,fYCOG,Geometry.ECal_front_z);
  double fE = sqrt(fMe*fMe + fBeamMomentum*fBeamMomentum);
  double fBG = fBeamMomentum/fSqrts; // beta gamma
  double fGam = sqrt(fBG*fBG+1.);
  double fBeta = fBG/fGam;

  //I think this is actually beta, not momentum
  //  TVector3 fBoostMom; fBoostMom.SetXYZ(fCOGAtECal.X()-fRTarg.X(),fCOGAtECal.Y()-fRTarg.Y(),fCOGAtECal.Z()-fRTarg.Z());
  TVector3 fBoostMom; fBoostMom.SetXYZ(fCOGAtECal.X()-fRTarg.X(),fCOGAtECal.Y()-fRTarg.Y(),deltaZ);
  fBoostMom *= (fBeta/fBoostMom.Mag());

  // if K = RMax/D is the max tangent in the lab, pi/2 - t < q*/2 < t, where t = atan(gamma RMax/D) = CoM angle must be > pi/4
  // t = pi/4 if gam = 1/K, i.e. at ~ 150 MeV

  // Maximum lab-frame angle
  double tanQMax = fRadiusMax/deltaZ;

  // Convert to center-of-mass frame angle using boost:
  // tan(theta*) ≈ gamma * tan(theta_lab)
  double tMax = TMath::ATan(fGam*tanQMax); //maximum angle in CoM frame

  //in CoM, particles have to be equal & opposite => angles add to pi/2 and neither can be > pi/4
  if (tMax < TMath::Pi()*0.25) {
    std::cout << "No solution? " << tMax << " " << fGam << " " << fBeamMomentum << std::endl;
  }

  //  double tanQMin = 1./(fGam*fGam*tanQMax);
  double tMin = 0.5*TMath::Pi() - tMax; //in CoM so angles add to pi/2
  
  return (deltaZ*TMath::Tan(tMin)/fGam);// (fCOGAtECal.Z()-fRTarg.Z())*tanQMin; 
}

double ECalResolution(double energy)
{
  //from ECal paper DOI 10.1088/1748-0221/15/10/T10003
  double firstterm = 2e-2/TMath::Sqrt(energy);
  double secondterm = 0.003e-2/energy;
  double thirdterm = 1.1e-2;
  
  double sigmaEoverE = TMath::Sqrt(firstterm*firstterm + secondterm*secondterm + thirdterm*thirdterm);
  double sigmeE = sigmaEoverE*energy;

  return sigmeE;
}

TVector3 VertexPosition(std::map<TString, TH2F*>& mh2, TString tag)
{
  double x = rand3.Gaus(0,Geometry.beamspotsize_X);
  double y = rand3.Gaus(0,Geometry.beamspotsize_Y);

  mh2.at("h" + tag + "VertexPosition")->Fill(x,y);
  
  return TVector3(x,y,Geometry.target_z);
}

//simulate tracker without B field, use energy from ECal
double RecoSqrtsStraightTracks(TLorentzVector pos_4mom, TLorentzVector ele_4mom, bool ECalSmearingOn)
{  
  //find transverse directions
  Double_t transverse_dir_pos_x = pos_4mom.X()/pos_4mom.Z();
  Double_t transverse_dir_pos_y = pos_4mom.Y()/pos_4mom.Z();
  Double_t transverse_dir_ele_x = ele_4mom.X()/ele_4mom.Z();
  Double_t transverse_dir_ele_y = ele_4mom.Y()/ele_4mom.Z();
      
  //positions at first scintillator positions
  Double_t x1_pos = transverse_dir_pos_x*(z1-Geometry.target_z);
  Double_t y1_pos = transverse_dir_pos_y*(z1-Geometry.target_z);
  Double_t x1_ele = transverse_dir_ele_x*(z1-Geometry.target_z);
  Double_t y1_ele = transverse_dir_ele_y*(z1-Geometry.target_z);

  Double_t x2_pos = transverse_dir_pos_x*(z2-Geometry.target_z);
  Double_t y2_pos = transverse_dir_pos_y*(z2-Geometry.target_z);
  Double_t x2_ele = transverse_dir_ele_x*(z2-Geometry.target_z);
  Double_t y2_ele = transverse_dir_ele_y*(z2-Geometry.target_z);
    
  //reconstruct momentum direction
  TVector3 reco_dir_pos = TVector3(x2_pos-x1_pos, y2_pos-y1_pos, z2-z1);
  reco_dir_pos = reco_dir_pos.Unit(); //normalise, this gives you only the direction
  TVector3 reco_dir_ele = TVector3(x2_ele-x1_ele, y2_ele-y1_ele, z2-z1);
  reco_dir_ele = reco_dir_ele.Unit();

  double Epos_true = pos_4mom.E();
  double Eele_true = ele_4mom.E();
  
  double Epos, Eele;
  if(!ECalSmearingOn)
    {
      Epos = Epos_true;
      Eele = Eele_true;
    }
  else
    {
      Double_t MaxECalRes_pos  = ECalResolution(Epos_true);
      Double_t MaxECalRes_ele  = ECalResolution(Eele_true);
      Double_t DeltaE_smear_pos = rand3.Uniform(-1.*MaxECalRes_pos,MaxECalRes_pos);
      Double_t DeltaE_smear_ele = rand3.Uniform(-1.*MaxECalRes_ele,MaxECalRes_ele);
      Epos = Epos_true+DeltaE_smear_pos;
      Eele = Eele_true+DeltaE_smear_ele;
    }
  
  double ptot_pos = TMath::Sqrt(Epos*Epos-me*me);
  double ptot_ele = TMath::Sqrt(Eele*Eele-me*me);
  
  //reconstruct 3-momenta
  TVector3 reco_3mom_pos = reco_dir_pos*ptot_pos;
  TVector3 reco_3mom_ele = reco_dir_ele*ptot_ele;
  
  //      std::cout<<"pos reco: "<<reco_3mom_pos_nores.Z()<<" "<<" pos nores: "<< pz_pos_lab <<" ele reco: "<<reco_3mom_ele_nores.Z()<<" ele nores: "<< pz_ele_lab<<std::endl;
  
  //reconstruct 4-momenta
  TLorentzVector reco_4mom_pos = TLorentzVector(reco_3mom_pos, Epos);
  TLorentzVector reco_4mom_ele = TLorentzVector(reco_3mom_ele, Eele);
  
  //total 4-momentum in lab frame = sum of e+/e- 4-momenta
  TLorentzVector reco_4mom_lab = reco_4mom_pos+reco_4mom_ele;

  return reco_4mom_lab.M();
}

double getBYfield(double x0, double y0, double z0in){ // input in mm
  //  double fConstantMagneticField = 0.;

  double fConstantMagneticFieldXmin = -26.0;//[cm]*cm;
  double fConstantMagneticFieldXmax =  26.0;//[cm]*cm;

  double fConstantMagneticFieldZmin = Geometry.magnet_front_z*0.1;//-37.5;// [cm]*cm;
  double fConstantMagneticFieldZmax = Geometry.magnet_back_z*0.1;//37.5;// [cm]*cm;

  double fSigmaFront = 27.4;//*cm;   //Based of the LNF magnetic measurement M. Raggi .ppt nov 2018
  double fSigmaBack  = 27.4;//*cm;   //Based of the LNF magnetic measurement M. Raggi .ppt nov 2018

  // The magnetic volume is a box centered at magnet center with x and y dimensions
  // as the magnet cavity and with z extends 50cm outside both sides of the magnet
  double fMagneticVolumePosZ = 0.;
  double fMagneticVolumeLengthX = 112.;//*cm;
  double fMagneticVolumeLengthY = 23.;//*cm;
  double fMagneticVolumeLengthZ = 100;//200.;//*cm;

  double x = x0*0.1;
  double y = y0*0.1;
  double z = z0in*0.1-fMagneticVolumePosZ;

  double BField0 = 1.;
  if ( (x<-0.5*fMagneticVolumeLengthX) || (x>0.5*fMagneticVolumeLengthX) ||
       (y<-0.5*fMagneticVolumeLengthY) || (y>0.5*fMagneticVolumeLengthY) ||
       (z<-0.5*fMagneticVolumeLengthZ) || (z>0.5*fMagneticVolumeLengthZ) ) {
    // Field outside magnetic volume is always null
    BField0 = 0.;
  } else if (x<fConstantMagneticFieldXmin || x>fConstantMagneticFieldXmax) {
    // Will need a function/map to smoothly send B0 to 0 along X
    BField0 = 0.;
  } else if (z<fConstantMagneticFieldZmin) {
    // Use gaussian to model upstream magnetic field rise
    double dZS = (z-fConstantMagneticFieldZmin)/fSigmaFront;
    BField0 = exp(-dZS*dZS);
  } else if (z>fConstantMagneticFieldZmax) {
    // Use gaussian to model downstream magnetic field fall
    double dZS = (z-fConstantMagneticFieldZmax)/fSigmaBack;
    BField0 = exp(-dZS*dZS);
  }
  return -fConstantMagneticField*BField0;
}

std::vector<TVector3> Swimmer(TVector3 LocalMomentum, TVector3 LocalPosition, std::map<TString, TH1F*>& mh1, std::map<TString, TH2F*>& mh2, const TString tag, const bool isOutsideDivergence = 0, const bool isOutsideBrem = 0)
{
  TVector3 PT(0,0,0);
  TVector3 LocalPositionRot(0,0,0);
  TVector3 magField(0,1,0); // T
  TVector3 magFieldS(0,1,0);// T m
	
  bool insideECal = kFALSE; // it is in the ECal
  int channelFired = -1;
  int ecalchannelFired = -1;
  double totalLength = 0;
  double ecalfireden = -1;
  double localpositionInECal[2] = {0,0};
  
  const double stepL = 0.2;//5;//mm 0.2;//0.1; // mm //without B field or fit, 10 mm for 1e5 particles = 4s, 1 mm for 1e5 particles = 28s, 0.1 mm for 1e5 particles = 280s
  double stepLM = stepL*1E-3; //m

  int ncells = Geometry.nECal_cells;
  
  const double bScale = 1;//0; // scale of the magnetic field

  std::vector<TVector3> vCoords; //vector to contain coordinates of partiles at tracker planes
  int charge = 0;
  if(tag == "Pos"||tag.Contains("Beam")||tag == "Brem") charge = +1;
  else if (tag == "Ele") charge = -1;
  else
    {
      std::cout<<"Swimmer: tag "<<tag<<" unclear"<<std::endl;
      return vCoords;
    }
  
  // swim loop
  
  int nSteps = (Geometry.ECal_front_z + Geometry.ECal_length - Geometry.target_z)/stepL*20;//100;
  for (int is = 0; is< nSteps; is++){
    magField.SetY(getBYfield(LocalPosition.X(),LocalPosition.Y(),LocalPosition.Z())*bScale);
    if (magField.Mag() > 1E-8){
      magFieldS.SetXYZ(magField.X()*stepLM,magField.Y()*stepLM,magField.Z()*stepLM); // Tm
      PT = LocalMomentum.Cross(magFieldS);
      PT *= charge*0.3*1E3/LocalMomentum.Mag(); // Length in m, PT in MeV //0.3 = conversion between SI (T, m) and natural units (GeV), 1e3 = GeV->MeV
      double localMag = LocalMomentum.Mag();
      LocalMomentum += PT;
      LocalMomentum *= (localMag/LocalMomentum.Mag());
    }
    TVector3 dStep(stepL*LocalMomentum.X()/LocalMomentum.Mag(),
		   stepL*LocalMomentum.Y()/LocalMomentum.Mag(),
		   stepL*LocalMomentum.Z()/LocalMomentum.Mag());
	  
    LocalPosition.SetX(LocalPosition.X() + dStep.X());
    LocalPosition.SetY(LocalPosition.Y() + dStep.Y());
    LocalPosition.SetZ(LocalPosition.Z() + dStep.Z());

    totalLength += stepL;
    
    // if (is % 500 == 0) {
	// std::cout << "z=" << LocalPosition.Z();
	// std::cout << " x=" << LocalPosition.X();
	// std::cout << " px=" << LocalMomentum.X();
	// std::cout << std::endl;
    //     }

    double prevZ = LocalPosition.Z() - dStep.Z();
    double x = LocalPosition.X();
    double y = LocalPosition.Y();

    //make plots at Z1
    if ( (prevZ < z1) && (LocalPosition.Z() >= z1) ) {
      FillPositionHistos(mh1, mh2, tag, x, y, "Z1", "");

      if(isOutsideDivergence) //for events where at least one particle has theta > beamdivergence
	FillPositionHistos(mh1, mh2, tag, x, y, "Z1", "BeamDivCut");
	
      if(isOutsideBrem) //for events where positron has theta > ThetaBremPos
	  FillPositionHistos(mh1, mh2, tag, x, y, "Z1", "BremPosCut");
      
      if(vCoords.size()!=0)
	{
	  std::cout<<"[Swimmer]: vCoords has non-zero size at Z1"<<std::endl;
	  break;
	}
      else vCoords.push_back(LocalPosition);
    }

    //make plots at Z2
    if ( (prevZ < z2) && (LocalPosition.Z() >= z2) ) {
      FillPositionHistos(mh1, mh2, tag, x, y, "Z2", "");

      if(isOutsideDivergence) //for events where at least one particle has theta > beamdivergence
	FillPositionHistos(mh1, mh2, tag, x, y, "Z2", "BeamDivCut");
	
      if(isOutsideBrem) //for events where positron has theta > ThetaBremPos
	  FillPositionHistos(mh1, mh2, tag, x, y, "Z2", "BremPosCut");

      if(vCoords.size()!=1)
	{
	  std::cout<<"[Swimmer]: vCoords has size != 1 at Z2"<<std::endl;
	  break;
	}
      else vCoords.push_back(LocalPosition);
    }
    
    if (LocalPosition.Z() > Geometry.ECal_front_z + Geometry.ECal_length) break;
    
    if (!insideECal) {
      if (LocalPosition.Z()-Geometry.ECal_front_z > 0 && LocalPosition.Z()-Geometry.ECal_front_z < Geometry.ECal_length) { // register position of entering in the ecal
	      
	int icellX = LocalPosition.X()/Geometry.ECal_cell_size +0.5 + ncells/2;
	int icellY = LocalPosition.Y()/Geometry.ECal_cell_size+0.5 + ncells/2;
	      
	if (LocalPosition.Perp() < Geometry.ECal_edge) {
		
	  if (icellX < ncells && !(icellX < 0)) { 
	    if (icellY < ncells && !(icellY < 0)) {
	      bool hole = TMath::Abs(icellX-ncells/2) <= 2 && TMath::Abs(icellY-ncells/2) <=2;
	      if (!hole) {
		insideECal = kTRUE;
		ecalchannelFired = icellX + ncells*icellY;
		ecalfireden = LocalMomentum.Mag();
		localpositionInECal[0] = LocalPosition.X();
		localpositionInECal[1] = LocalPosition.Y();
		break;
	      }
	    }
	  }
	}
      }
    }
  }
  return vCoords;
}

void BeamSpotSpread(TString tag, std::map<TString, TH1F*>& mh1, std::map<TString, TH2F*>& mh2)
{
  double sqrts_sideband = 0;
  if(tag == "Beam_LowSideband") sqrts_sideband = 16.5;
  else if(tag == "Beam_HighSideband") sqrts_sideband = 20.0;
  else
    {
      std::cout<<"BeamSpotSpread: tag "<<tag<<" unknown"<<std::endl;
      return;
    }  
  int nparticles = 1e3;
  double Ebeam_sideband = (sqrts_sideband*sqrts_sideband)/(2*me)-me;
  double pbeam_sideband = TMath::Sqrt(Ebeam_sideband*Ebeam_sideband-me*me);

  if(tag.Contains("Low"))
    mh1["hPbeamSmear_LowSideband"] = new TH1F("hPbeamSmear_LowSideband","hPbeamSmear_LowSideband",100,pbeam_sideband-5,pbeam_sideband+5);
  else if(tag.Contains("High"))
    mh1["hPbeamSmear_HighSideband"] = new TH1F("hPbeamSmear_HighSideband","hPbeamSmear_HighSideband",100,pbeam_sideband-5,pbeam_sideband+5);
  
  std::cout<<pbeam_sideband<<std::endl;
  
  for (int ii = 0; ii<nparticles; ii++)
    {
      //smear beam energy
      Double_t SigmaPbeam = SigmaPOverP*pbeam_sideband;
      Double_t DeltaP = rand3.Gaus(0, SigmaPbeam);
      Double_t PbeamSmear = pbeam_sideband+DeltaP;
      if(tag.Contains("Low"))
	mh1.at("hPbeamSmear_LowSideband")->Fill(PbeamSmear);
      else if(tag.Contains("High"))
	mh1.at("hPbeamSmear_HighSideband")->Fill(PbeamSmear);
      
      Swimmer(TVector3(0,0,PbeamSmear), VertexPosition(mh2, "Beam"), mh1, mh2, tag);
    }
}

void DrawAllHists(std::map<TString, TH1F*>& mh1, std::map<TString, TH2F*>& mh2)
{
  //  TFile* fOut = new TFile(sOutputFilePath, "recreate");

  // 1D histograms
  for (auto& [name, hist] : mh1) {
    TCanvas* c = new TCanvas("c_" + name, name, 900, 700);
    hist->Draw();
    //    hist->Write();
  }
  
  // 2D histograms
  for (auto& [name, hist] : mh2) {
    TCanvas* c = new TCanvas("c_" + name, name, 900, 700);
    hist->Draw("colz");
    //    hist->Write();
  }
  // fOut->Write();
  // fOut->Close();
}

//reconstructed momentum from tracks bent in B field
double RecoMomCurvedTrack(TVector3 Coord_Z1, TVector3 Coord_Z2)
{
  //p[GeV] = 0.3qBL/theta

  double X1 = Coord_Z1.X(); double X2 = Coord_Z2.X(); double DeltaX = X2-X1;
  double Y1 = Coord_Z1.Y(); double Y2 = Coord_Z2.Y(); double DeltaY = Y2-Y1;
  double Z1 = Coord_Z1.Z(); double Z2 = Coord_Z2.Z(); double DeltaZ = Z2-Z1;
  
  double tantheta = TMath::Sqrt(DeltaX*DeltaX+DeltaY*DeltaY)/DeltaZ;
  double theta = TMath::ATan(tantheta);

  double mag_field_length = Geometry.magnet_back_z-Geometry.magnet_front_z;
  
  double momentum = 0.3*fConstantMagneticField*mag_field_length/theta;
  return momentum;
}

double RecoSqrtsCurvedTracks(std::vector<TVector3> vPosCoord, std::vector<TVector3> vEleCoord)
{
  TVector3 pos_atZ1 = vPosCoord[0]; TVector3 pos_atZ2 = vPosCoord[1];
  TVector3 ele_atZ1 = vEleCoord[0]; TVector3 ele_atZ2 = vEleCoord[1];
  
  std::cout<<RecoMomCurvedTrack(pos_atZ1, pos_atZ2)<<" "<<RecoMomCurvedTrack(ele_atZ1, ele_atZ2)<<std::endl;
  double sqrts;
  return sqrts;
}

double TrackerResolution(double momentum) //momentum in MeV
{
  double GeVMomentum = momentum*1e-3;

  double mag_field_length = 0.75;//(Geometry.magnet_back_z-Geometry.magnet_front_z)*1e-3; //L, m
  double plane_separation_length = 0.5;//(Geometry.chamber_back_z-Geometry.magnet_front_z)*1e-3/nTrackerPlanes; //l, m

  double nTrackerDetectors = 2*nTrackerPlanes;
  
  double relative_resolution = 8/0.3*1/(fConstantMagneticField*mag_field_length)*position_res/(mag_field_length*TMath::Sqrt(nTrackerDetectors))*GeVMomentum;
  double momentum_resolution = relative_resolution*momentum; //MeV

  //  std::cout<<"mag_field_length "<<mag_field_length<<" plane_separation_length "<<plane_separation_length<<" nTrackerPlanes "<<nTrackerPlanes<<" fConstantMagneticField "<<fConstantMagneticField<<" GeVMomentum "<<GeVMomentum<<" resolution "<<relative_resolution<<std::endl;
  
  //  std::cout<<"momentum "<<momentum<<" relative_resolution "<<relative_resolution<<" momentum_resolution "<<momentum_resolution<<std::endl;
  return momentum_resolution;
}

TVector3 MomentumBuilder(double momentum, double cos_theta, double phi)
{
  //trigonometry
  Double_t theta = TMath::ACos(cos_theta);
  Double_t sin_theta = TMath::Sin(theta);
  Double_t sin_phi   = TMath::Sin(phi);
  Double_t cos_phi   = TMath::Cos(phi);
  
  //momentum components in CoM
  Double_t pz = momentum*cos_theta;
  Double_t px = momentum*sin_theta*cos_phi;
  Double_t py = momentum*sin_theta*sin_phi;

  return TVector3(px, py, pz);
}

void InitialiseHistos(std::map<TString, TH1F*>& mh1, std::map<TString, TH2F*>& mh2)
{
  mh1["hCosThetaPosCoMPassing"] = new TH1F("hCosThetaPosCoMPassing",";hCosThetaPosCoMPassing;",100,-1,1);

  mh1["hPos_mom"] = new TH1F("hPos_mom","momentum of positron from decay", 300, 0, 300);
  mh1["hEle_mom"] = new TH1F("hEle_mom","momentum of electron from decay", 300, 0, 300);
  
  double xlim = Geometry.magnet_width_x/2.;
  double ylim = Geometry.magnet_width_y/2.;

  for(const TString& tag : aTags)
    {
      TString sTagTitle;
      if(tag == "Beam") sTagTitle = "beam";
      else if(tag == "Brem") sTagTitle = "bremsstrahlung";
      else if(tag == "Ele"||tag == "Pos") sTagTitle = "decay";
      else sTagTitle = "unkown tag";
      
      for(const auto& [sZKey, zPos] : vpZPositions)
	{
	  TString particle = "electron";
	  if(tag == "Pos") particle = "positron";

	  //initalise 2D position histograms
	  //basic 2D histograms before any acceptance cuts
	  TString base2dhistoname = Form("h%s_XY%s", tag.Data(), sZKey.Data());	      
	  TString base2dhistotitle = Form("position of %s from %s at %s", particle.Data(), sTagTitle.Data(), sZKey.Data());
	  mh2[base2dhistoname] = new TH2F(base2dhistoname, base2dhistotitle, 200, -xlim, xlim, 200, -ylim, ylim);

	  //	  if(tag!="Ele" && tag!="Pos") continue; //cuts are only interesting for decays
	  
	  //2D histograms including acceptance cuts
	  for(const auto& [sCutName, sCutTitle] : vpsCuts)
	    {
	      TString cut2dhistoname = Form("h%s_XY%s_%s", tag.Data(), sZKey.Data(), sCutName.Data());
	      TString cut2dhistotitle = Form("position of %s from %s at %s %s", particle.Data(), sTagTitle.Data(), sZKey.Data(), sCutTitle.Data());
	      mh2[cut2dhistoname] = new TH2F(cut2dhistoname, cut2dhistotitle, 200, -xlim, xlim, 200, -ylim, ylim);
	    }
	  
	  //initalise 1D position histograms
	  for(const TString& sXY : aXY)
	    {
	      //basic 1D histograms before any acceptance cuts
	      TString base1dhistoname = Form("h%s_%s%s", tag.Data(), sXY.Data(), sZKey.Data());	      
	      TString base1dhistotitle = Form("%s position of %s from decay at %s", sXY.Data(), particle.Data(), sZKey.Data());
	      mh1[base1dhistoname] = new TH1F(base1dhistoname, base1dhistotitle, 200, -xlim, xlim);
	      
	      //1D histograms including acceptance cuts
	      for(const auto& [sCutName, sCutTitle] : vpsCuts)
		{
		  TString cut1dhistoname = Form("h%s_%s%s_%s", tag.Data(), sXY.Data(), sZKey.Data(), sCutName.Data());
		  TString cut1dhistotitle = Form("%s position of %s from decay at %s %s", sXY.Data(), particle.Data(), sZKey.Data(), sCutTitle.Data());
		  
		  mh1[cut1dhistoname] = new TH1F(cut1dhistoname, cut1dhistotitle, 200, -xlim, xlim);
		}
	    }
       	}
    }

  mh1["hRelativeMomResTrack"] = new TH1F("hRelativeMomResTrack","",100,0,0.1);

  mh1["hSqrtsTrue"]           = new TH1F("hSqrtsTrue","",100,16,18);
  mh1["hSqrtsTrackerRes"]     = new TH1F("hSqrtsTrackerRes",";reconstructed #sqrt{s}",200,16,18);
  mh1["hDeltaSqrtsNoRes"]     = new TH1F("hDeltaSqrtsNoRes",";#Delta(#sqrt{s})_{True} (MeV)",100,-2,2);
  mh1["hDeltaSqrtsECalSmear"] = new TH1F("hDeltaSqrtsECalSmear",";#Delta(#sqrt{s})_{ECalSmear} (MeV)",100,-2,2);

  // 2D histograms
  mh2["hBeam_LowSideband_XYZ1"]  = new TH2F("hBeam_LowSideband_XYZ1","Position of uninteracted beam at low sideband energy at z1",50,-1,1,50,-1,1);
  mh2["hBeam_LowSideband_XYZ2"]  = new TH2F("hBeam_LowSideband_XYZ2","Position of uninteracted beam at low sideband energy at z2",100,128,132,10,-1,1);
  mh2["hBeam_HighSideband_XYZ1"] = new TH2F("hBeam_HighSideband_XYZ1","Position of uninteracted beam at high sideband energy at z1",50,-1,1,50,-1,1);
  mh2["hBeam_HighSideband_XYZ2"] = new TH2F("hBeam_HighSideband_XYZ2","Position of uninteracted beam at high sideband energy at z2",100,86,90,100,-1,1);

  mh2["hDecayVertexPosition"] = new TH2F("hDecayVertexPosition","Position of decay vertex on target",200,-20,20,200,-20,20);
  mh2["hBremVertexPosition"] = new TH2F("hBremVertexPosition","Position of bremsstrahlung vertex on target",200,-20,20,200,-20,20);

  mh1["hBeam_LowSideband_XZ1"]  = new TH1F("hBeam_LowSideband_XZ1","X position of uninteracted beam at low sideband energy at z1",100,-1,1);
  mh1["hBeam_LowSideband_XZ2"]  = new TH1F("hBeam_LowSideband_XZ2","X position of uninteracted beam at low sideband energy at z2",100,125,135);
  mh1["hBeam_HighSideband_XZ1"] = new TH1F("hBeam_HighSideband_XZ1","X position of uninteracted beam at high sideband energy at z1",100,-1,1);
  mh1["hBeam_HighSideband_XZ2"] = new TH1F("hBeam_HighSideband_XZ2","X position of uninteracted beam at high sideband energy at z2",100,86,90);

  mh1["hBeam_LowSideband_YZ1"]  = new TH1F("hBeam_LowSideband_YZ1","Y position of uninteracted beam at low sideband energy at z1",100,-1,1);
  mh1["hBeam_LowSideband_YZ2"]  = new TH1F("hBeam_LowSideband_YZ2","Y position of uninteracted beam at low sideband energy at z2",100,-1,1);
  mh1["hBeam_HighSideband_YZ1"] = new TH1F("hBeam_HighSideband_YZ1","Y position of uninteracted beam at high sideband energy at z1",100,-1,1);
  mh1["hBeam_HighSideband_YZ2"] = new TH1F("hBeam_HighSideband_YZ2","Y position of uninteracted beam at high sideband energy at z2",100,-1,1);

  mh2["hThetaVsEnergy"] = new TH2F("hThetaVsEnergy","hThetaVsEnergy",200,0,0.2,250,0,500);
}

void TrackerToyMC()
{
  TBenchmark Benchmark = TBenchmark();
  Benchmark.Start("macro");

  //set up histos
  std::map<TString, TH1F*> mh1;
  std::map<TString, TH2F*> mh2;
  InitialiseHistos(mh1, mh2);

  TString sBremBinFile = "/Users/elizabeth/Documents/PADME/PADME-Tracker-Toy/BremSampler.bin";
  
  int nDecays = 1e5;

  int nOutsideDivergence = 0; //no. decays where both decay products have theta > beam divergence
  int nOutsideBrem       = 0; //no. decays where decay positron has theta > ThetaBremPos

  Double_t PbeamNom = 279; //nominal beam mom, before adding energy spread, in MeV (scan between 262 MeV and 296 MeV gives 279 MeV average)
  
  for(int ii = 0; ii<nDecays; ii++)
    {      
      bool isOutsideDivergence = 0; //both decay products have theta > beam divergence
      bool isOutsideBrem = 0;       //decay positron has theta > ThetaBremPos

      double Pbeam = PbeamNom;
      if(ii%1000 == 0) std::cout<<"particle "<<ii<<std::endl;
      if(BeamEnergySpread) Pbeam = Pbeam+rand3.Gaus(0,SigmaPOverP*Pbeam);
      
      //set beam conditions
      Double_t Ebeam = TMath::Sqrt(Pbeam*Pbeam+me*me); // MeV
      Double_t sqrts = TMath::Sqrt(2*me*(Ebeam+me));
      
      double signalweight = 1;
      //      if(ElectronEnergySpread) signalweight = SignalShape(13.6,0.0025,1.7,sqrts,sqrts); //from conversation with Tommaso
      
      TLorentzVector system(0,0,Pbeam,Ebeam+me); //first approximation - electron energy in lab frame is negligible wrt beam energy
      TVector3 boost = system.BoostVector();
      
      //set acceptance limits of ECal
      Double_t RMax = Geometry.ECal_RMax;                        //max radius of highest energy ECal cluster
      /*
      Double_t RMin = kineRMin(sqrts);                           //min distance of cluster
      Double_t tanmin = RMin/Geometry.deltaZ;
      Double_t tanmax = RMax/Geometry.deltaZ;
      */
      //generate decay in CoM frame
      Double_t Ee_CoM = 0.5*sqrts;                      //Energy of e+/e- in CoM is half of sqrts (equal and opposite)
      Double_t mom_CoM = TMath::Sqrt(Ee_CoM*Ee_CoM-me*me); //magnitude of momentum of e+/e- in CoM
      
      //theta = polar angle: angle wrt z axis, phi = azimuthal angle: rotation around z axis
      Double_t cos_theta_pos_CoM = rand3.Uniform(-1,1);   //polar angle in CoM is uniformly distributed in cos(theta)
      Double_t phi_pos = rand3.Uniform(0,2.*TMath::Pi()); //azimuthal angle in CoM is uniformly distributed between 0-2pi

      TVector3 pos_3mom_CoM = MomentumBuilder(mom_CoM, cos_theta_pos_CoM, phi_pos);
      TVector3 ele_3mom_CoM = -1.*pos_3mom_CoM; //3-momentum of positron in CoM
      
      TLorentzVector pos_4mom_CoM = TLorentzVector(pos_3mom_CoM,Ee_CoM);  //4-momentum of positron in CoM
      TLorentzVector ele_4mom_CoM = TLorentzVector(-pos_3mom_CoM,Ee_CoM); //4-momentum of electron in CoM is equal and opposite to positron
      
      TLorentzVector pos_4mom_lab = pos_4mom_CoM;
      pos_4mom_lab.Boost(boost);
      TLorentzVector ele_4mom_lab = ele_4mom_CoM;
      ele_4mom_lab.Boost(boost);

      //      std::cout<<"theta pos "<<pos_4mom_lab.Theta()*1e3<<" ele "<<ele_4mom_lab.Theta()<<std::endl;
      
      if(pos_4mom_lab.Theta()*1e3 > Geometry.beamdivergence && ele_4mom_lab.Theta()*1e3 > Geometry.beamdivergence)
	{
	  isOutsideDivergence = 1;
	  nOutsideDivergence++;
	}
      
      if(pos_4mom_lab.Theta()*1e3 > ThetaBremPos)
	{
	  isOutsideBrem = 1;
	  nOutsideBrem++;
	}

      //total energy and mom of e+/e- in lab
      Double_t ptot_pos_lab = pos_4mom_lab.P();
      Double_t ptot_ele_lab = ele_4mom_lab.P();
      Double_t Epos_true = pos_4mom_lab.E();
      Double_t Eele_true = ele_4mom_lab.E();

      mh1.at("hPos_mom")->Fill(ptot_pos_lab);
      mh1.at("hEle_mom")->Fill(ptot_ele_lab);
      
      //find momentum resolution from tracker: 
      double pos_track_res = TrackerResolution(ptot_pos_lab);
      double ele_track_res = TrackerResolution(ptot_ele_lab);
      mh1["hRelativeMomResTrack"]->Fill(pos_track_res/ptot_pos_lab);
      mh1["hRelativeMomResTrack"]->Fill(ele_track_res/ptot_ele_lab);
	
      //for straight track approack, check for acceptance as soon as you've boosted
      Double_t pt_pos_lab = pos_4mom_lab.Pt();
      Double_t pt_ele_lab = ele_4mom_lab.Pt();
      Double_t tan_theta_pos_lab = pt_pos_lab/pos_4mom_lab.Z();
      Double_t tan_theta_ele_lab = pt_ele_lab/ele_4mom_lab.Z();
      
      //      std::cout<<tan_theta_pos_lab<<" "<<tan_theta_ele_lab<<" "<<tanmin<<" "<<tanmax<<std::endl;

      double reco_sqrts_straight_nores = -999;
      double reco_sqrts_straight_ECalSmear = -999;
      /*      if(!(tan_theta_pos_lab<tanmin||tan_theta_ele_lab<tanmin||tan_theta_pos_lab>tanmax||tan_theta_ele_lab>tanmax))
	{
	  mh1.at("hCosThetaPosCoMPassing")->Fill(cos_theta_pos_CoM);

	  reco_sqrts_straight_nores     = RecoSqrtsStraightTracks(pos_4mom_lab, ele_4mom_lab, 0);
	  mh1.at("hDeltaSqrtsNoRes")->Fill(reco_sqrts_straight_nores-sqrts);

	  reco_sqrts_straight_ECalSmear = RecoSqrtsStraightTracks(pos_4mom_lab, ele_4mom_lab, 1);
	  mh1.at("hDeltaSqrtsECalSmear")->Fill(reco_sqrts_straight_ECalSmear-sqrts);
	}
      */
      double reco_sqrts_straight_truth     = RecoSqrtsStraightTracks(pos_4mom_lab, ele_4mom_lab, 0);
      mh1["hSqrtsTrue"]->Fill(reco_sqrts_straight_truth,signalweight);
      
      TVector3 vertex_position = VertexPosition(mh2, "Decay");
      TVector3 pos_3mom_lab_true = TVector3(pos_4mom_lab.X(), pos_4mom_lab.Y(), pos_4mom_lab.Z());
      TVector3 ele_3mom_lab_true = TVector3(ele_4mom_lab.X(), ele_4mom_lab.Y(), ele_4mom_lab.Z());

      //smearing from tracker
      double ptot_pos_lab_smear = ptot_pos_lab + rand3.Gaus(0,pos_track_res);
      double ptot_ele_lab_smear = ptot_ele_lab + rand3.Gaus(0,ele_track_res);
      
      TVector3 dir_pos_true = pos_3mom_lab_true.Unit();
      TVector3 dir_ele_true = ele_3mom_lab_true.Unit();
      
      TVector3 reco_3mom_pos_smear = dir_pos_true * ptot_pos_lab_smear;
      TVector3 reco_3mom_ele_smear = dir_ele_true * ptot_ele_lab_smear;
      
      double Epos_smear = TMath::Sqrt(ptot_pos_lab_smear*ptot_pos_lab_smear + me*me);
      double Eele_smear = TMath::Sqrt(ptot_ele_lab_smear*ptot_ele_lab_smear + me*me);
      
      TLorentzVector reco_4mom_pos_smear = TLorentzVector(reco_3mom_pos_smear,Epos_smear);
      TLorentzVector reco_4mom_ele_smear = TLorentzVector(reco_3mom_ele_smear,Eele_smear);
      
      TLorentzVector reco_4mom_lab_smear = reco_4mom_pos_smear+reco_4mom_ele_smear;
      
      double reco_sqrts_track_smear = reco_4mom_lab_smear.M();
      mh1["hSqrtsTrackerRes"]->Fill(reco_sqrts_track_smear);
      
      //swimmer
      std::vector<TVector3> vPosCoord = Swimmer(pos_3mom_lab_true, vertex_position, mh1, mh2, "Pos", isOutsideDivergence, isOutsideBrem);
      std::vector<TVector3> vEleCoord = Swimmer(ele_3mom_lab_true, vertex_position, mh1, mh2, "Ele", isOutsideDivergence, isOutsideBrem);

      if(vPosCoord.size()!= vpZPositions.size() || vEleCoord.size()!= vpZPositions.size())
	{
	  std::cout<<"ii "<<ii<<" vPosCoord has size "<<vPosCoord.size()<<" "<<" vEleCoord has size "<<vEleCoord.size()<<" "<<std::endl;
	  continue;
	}

      double reco_sqrts_curved = -999;
      //      RecoSqrtsCurvedTracks(vPosCoord, vEleCoord);
    }//end loop over decays

  double nBrem = 1e4;
  for(int jj = 0; jj<nBrem; jj++)
    {
      if(jj%1000 == 0) std::cout<<"Bremsstrahlung "<<jj<<std::endl;
      std::pair<double, double> pThetaEne = BremThetaEneSampler(sBremBinFile, mh2);
      Double_t phi_brem = rand3.Uniform(0,2.*TMath::Pi()); //azimuthal angle in CoM is uniformly distributed between 0-2pi

      double theta  = pThetaEne.first;
      double energy = pThetaEne.second;
      
      double momentum = energy*energy-me*me;
      TVector3 bremMom = MomentumBuilder(momentum, TMath::Cos(theta), phi_brem);

      TVector3 vertex_position = VertexPosition(mh2, "Brem");
      
      std::vector<TVector3> vBremCoord = Swimmer(bremMom, vertex_position, mh1, mh2, "Brem");
    }
  

  // BeamSpotSpread("Beam_LowSideband", mh1, mh2);
  // BeamSpotSpread("Beam_HighSideband", mh1, mh2);
  
  Benchmark.Show("macro");

  DrawAllHists(mh1, mh2);

  int nelectronsim = 1e4;
  double sqrtsmin = 16;
  double sqrtsmax = 18;
  double sqrtsstep = (sqrtsmax-sqrtsmin)/nelectronsim;

  TGraph* gr = new TGraph(nelectronsim);
  
  for(int ii = 0; ii<nelectronsim; ii++)
    {
      double sqrts = sqrtsmin+ii*sqrtsstep;
      double sigheight = SignalShape(13.6,0.0025,1.7,16.9015,sqrts); //from conversation with Tommaso
      gr->SetPoint(ii, sqrts, sigheight);
    }


  TCanvas* cSqrts = new TCanvas("cSqrtsElectron","cSqrtsElectron",900,700);
  gr->SetTitle(";#sqrt{s} (MeV);Signal (aribtrary units)");
  gr->SetLineWidth(2);
  gr->Draw("AL");  // A = axes, L = line
  TH1* hHisto = mh1.at("hSqrtsTrackerRes");
  hHisto->Scale(30./nDecays);
  hHisto->Draw("histo same");

  // Create a custom TPaveText (no dividing lines)
  TPaveText *pt = new TPaveText(0.57,0.78,0.899,0.899,"NDC");
  pt->SetFillColor(0);          // transparent
  pt->SetLineColor(kRed);       // red border
  pt->SetLineWidth(3);          // border width
  pt->SetTextFont(42);          // clean font
  pt->SetTextSize(0.04);        // optional text size
  pt->SetTextAlign(32);         // right-aligned
  pt->SetBorderSize(1);         // width - no shadow

  // Add mean and sigma
  pt->AddText(Form("Electron Mean = %.3f", gr->GetMean()));
  pt->AddText(Form("Electron RMS = %.3f", gr->GetRMS()));

  // Draw the box
  pt->Draw();

  // --- Blue box for histogram ---
  TPaveText *ptBlue = new TPaveText(0.57,0.68,0.899,0.77,"NDC"); // slightly below red box
  ptBlue->SetFillColor(0);
  ptBlue->SetLineColor(kBlue);
  ptBlue->SetLineWidth(3);
  ptBlue->SetTextFont(42);
  ptBlue->SetTextSize(0.04);
  ptBlue->SetTextAlign(32);
  ptBlue->SetBorderSize(1);
  
  ptBlue->AddText(Form("Tracker Mean = %.3f", hHisto->GetMean()));
  ptBlue->AddText(Form("Tracker Std Dev = %.3f", hHisto->GetStdDev()));
  
  ptBlue->Draw();
  
  gPad->Modified();
  gPad->Update();


  std::cout<<"beam energy "<<PbeamNom<<" divergence accetpance "<<1.*nOutsideDivergence/nDecays<<" Brem acceptance "<<1.*nOutsideBrem/nDecays<<std::endl;
  
}
