struct PADMEGeometry
{
  //coordinates are wrt the centre of the PADME dipole
  double target_z       = -732; //mm, Run 4 position (run 3 = -1028 mm, run 4 = 510 mm)
  double magnet_front_z = -490; //mm
  double magnet_back_z  = 510;  //mm
  double chamber_back_z = 2300; //mm
  double ECal_front_z   = 2436; //mm
  double ECal_length    = 230;  //mm
  double magnet_width_x = 435;  //mm
  double magnet_width_y = 205;  //mm
  double ECal_RMax      = 270;  //mm, maximum radius of highest energy cluster
  double ECal_cell_size = 21;   //mm
  int    nECal_cells    = 29;   // per row or column

  double ECal_edge      = (nECal_cells/2+0.5)*ECal_cell_size; // 304.5mm
  double deltaZ = ECal_front_z-target_z; //target-ECal distance
};

PADMEGeometry Geometry;
const Double_t me = 0.511; //MeV
const Double_t z1 = Geometry.magnet_front_z; //first scintillator position at entrance to magnet
const Double_t z2 = Geometry.chamber_back_z; //second scintillator position at back of chamber
 
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
      Epos = Epos_true+ECalResolution(Epos_true);
      Eele = Eele_true+ECalResolution(Eele_true);
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
  double fConstantMagneticField = 0.36;//0.4542; // [T] 
  //  fConstantMagneticField = 0.;                                                                                                                                                                      

  double fConstantMagneticFieldXmin = -26.0;//[cm]*cm;
  double fConstantMagneticFieldXmax =  26.0;//[cm]*cm;

  double fConstantMagneticFieldZmin = -37.5;// [cm]*cm;
  double fConstantMagneticFieldZmax =  37.5;// [cm]*cm;

  double fSigmaFront = 27.4;//*cm;   //Based of the LNF magnetic measurement M. Raggi .ppt nov 2018                                                                                                             
  double fSigmaBack  = 27.4;//*cm;   //Based of the LNF magnetic measurement M. Raggi .ppt nov 2018                                                                                                             

  // The magnetic volume is a box centered at magnet center with x and y dimensions                                                                                                                   
  // as the magnet cavity and with z extends 50cm outside both sides of the magnet                                                                                                                    
  double fMagneticVolumePosZ = 0.;
  double fMagneticVolumeLengthX = 112.;//*cm;
  double fMagneticVolumeLengthY = 23.;//*cm;
  double fMagneticVolumeLengthZ = 200.;//*cm;

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

void Swimmer(TVector3 LocalDirection,  TVector3 LocalPosition, int charge)
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
  
  const double stepL = 10;//mm 0.2;//0.1; // mm
  double stepLM = stepL*1E-3; //m

  int ncells = Geometry.nECal_cells;
  
  const double bScale = 0; // scale of the magnetic field (Beth: T I guess?)
  
  // swim loop
  
  int nSteps = (Geometry.ECal_front_z + Geometry.ECal_length - Geometry.target_z)/stepL*20;//100;
  for (int is = 0; is< nSteps; is++){
    magField.SetY(getBYfield(LocalPosition.X(),LocalPosition.Y(),LocalPosition.Z())*bScale);
    if (magField.Mag() > 1E-8){
      magFieldS.SetXYZ(magField.X()*stepLM,magField.Y()*stepLM,magField.Z()*stepLM); // Tm
      PT = LocalDirection.Cross(magFieldS);
      PT *= charge*0.3*1E3/LocalDirection.Mag(); // Length in m, PT in MeV
      double localMag = LocalDirection.Mag();
      LocalDirection += PT;
      LocalDirection *= (localMag/LocalDirection.Mag());
    }
    TVector3 dStep(stepL*LocalDirection.X()/LocalDirection.Mag(),
		   stepL*LocalDirection.Y()/LocalDirection.Mag(),
		   stepL*LocalDirection.Z()/LocalDirection.Mag());
	  
    LocalPosition.SetX(LocalPosition.X() + dStep.X());
    LocalPosition.SetY(LocalPosition.Y() + dStep.Y());
    LocalPosition.SetZ(LocalPosition.Z() + dStep.Z());

    totalLength += stepL;

    if (/*fVerbose && */(is%30 == 0)) {
      cout << "In detailed magnetic Transport: PT kick = ( "<< PT.X() << ", "<< PT.Y()<<", "<< PT.Z()<< ") "
	   << " mom[MeV] = ( " << LocalDirection.X() << " ,"<<LocalDirection.Y() << " ,"<< LocalDirection.Z() << ")"
	   << " pos[mm] = ( " << LocalPosition.X() << " ,"<<LocalPosition.Y() << " ,"<< LocalPosition.Z() << ")"
	   << " local pos[mm] = ( " << LocalPositionRot.X() << " ,"<<LocalPositionRot.Y() << " ,"<< LocalPositionRot.Z() << ")"
	   << " magField[Tm] = ( " << magFieldS.X() << " ,"<<magFieldS.Y() << " ,"<<magFieldS.Z() << ")"<< endl;
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
		ecalfireden = LocalDirection.Mag();
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
}

void TrackerToyMC()
{
  //set beam conditions
  Double_t Pbeam = 279; //MeV (scan between 262 MeV and 296 MeV gives 279 MeV average)
  Double_t Ebeam = TMath::Sqrt(Pbeam*Pbeam+me*me); // MeV
  Double_t sqrts = TMath::Sqrt(2*me*(Ebeam+me));
  TLorentzVector system(0,0,Pbeam,Ebeam+me); //first approximation - electron energy in lab frame is negligible wrt beam energy
  TVector3 boost = system.BoostVector();

  //set acceptance limits of ECal
  Double_t RMax = Geometry.ECal_RMax;                        //max radius of highest energy ECal cluster
  Double_t RMin = kineRMin(sqrts);                           //min distance of cluster
  Double_t tanmin = RMin/Geometry.deltaZ;
  Double_t tanmax = RMax/Geometry.deltaZ;

  //set up histos
  TH1F* hCosThetaPosCoMPassing = new TH1F("hCosThetaPosCoMPassing",";hCosThetaPosCoMPassing;",100,-1,1);

  TH1F* hDeltaSqrtsNoRes       = new TH1F("hDeltaSqrtsNoRes"      ,";#Delta(#sqrt{s})_{True} (MeV)",100,-2,2);
  TH1F* hDeltaSqrtsECalSmear   = new TH1F("hDeltaSqrtsECalSmear"  ,";#Delta(#sqrt{s})_{ECalSmear} (MeV)",100,-2,2);
  
  int nparticles = 1e5;
  TRandom3 rand = TRandom3(1);
  
  for(int ii = 0; ii<nparticles; ii++)
    {
      //generate decay in CoM frame
      Double_t Ee_CoM = 0.5*sqrts;                      //Energy of e+/e- in CoM is half of sqrts (equal and opposite)
      Double_t mom_CoM = TMath::Sqrt(Ee_CoM*Ee_CoM-me*me); //magnitude of momentum of e+/e- in CoM
      
      //theta = polar angle: angle wrt z axis, phi = azimuthal angle: rotation around z axis
      Double_t cos_theta_pos_CoM = rand.Uniform(-1,1);   //polar angle in CoM is uniformly distributed in cos(theta)
      Double_t phi_pos = rand.Uniform(0,2.*TMath::Pi()); //azimuthal angle in CoM is uniformly distributed between 0-2pi

      //trigonometry
      Double_t theta_pos_CoM = TMath::ACos(cos_theta_pos_CoM);
      Double_t sin_theta_pos_CoM = TMath::Sin(theta_pos_CoM);
      Double_t sin_phi_pos_CoM   = TMath::Sin(phi_pos);
      Double_t cos_phi_pos_CoM   = TMath::Cos(phi_pos);

      //momentum components in CoM
      Double_t pz_pos_CoM = mom_CoM*cos_theta_pos_CoM;
      Double_t px_pos_CoM = mom_CoM*sin_theta_pos_CoM*cos_phi_pos_CoM;
      Double_t py_pos_CoM = mom_CoM*sin_theta_pos_CoM*sin_phi_pos_CoM;

      TVector3 pos_3mom_CoM = TVector3(px_pos_CoM, py_pos_CoM, pz_pos_CoM); //3-momentum of positron in CoM
      TVector3 ele_3mom_CoM = -1.*pos_3mom_CoM; //3-momentum of positron in CoM
      
      TLorentzVector pos_4mom_CoM = TLorentzVector(pos_3mom_CoM,Ee_CoM);  //4-momentum of positron in CoM
      TLorentzVector ele_4mom_CoM = TLorentzVector(-pos_3mom_CoM,Ee_CoM); //4-momentum of electron in CoM is equal and opposite to positron
      
      TLorentzVector pos_4mom_lab = pos_4mom_CoM;
      pos_4mom_lab.Boost(boost);
      TLorentzVector ele_4mom_lab = ele_4mom_CoM;
      ele_4mom_lab.Boost(boost);

      //total energy and mom of e+/e- in lab
      Double_t ptot_pos_lab = pos_4mom_lab.P();
      Double_t ptot_ele_lab = ele_4mom_lab.P();
      Double_t Epos_true = pos_4mom_lab.E();
      Double_t Eele_true = ele_4mom_lab.E();
      
      //for straight track approack, check for acceptance as soon as you've boosted
      Double_t pt_pos_lab = pos_4mom_lab.Pt();
      Double_t pt_ele_lab = ele_4mom_lab.Pt();
      Double_t tan_theta_pos_lab = pt_pos_lab/pos_4mom_lab.Z();
      Double_t tan_theta_ele_lab = pt_ele_lab/ele_4mom_lab.Z();

      std::cout<<tan_theta_pos_lab<<" "<<tan_theta_ele_lab<<" "<<tanmin<<" "<<tanmax<<std::endl;
      
      double reco_sqrts_straight_nores = -999;
      double reco_sqrts_straight_ECalSmear = -999;
      if(!(tan_theta_pos_lab<tanmin||tan_theta_ele_lab<tanmin||tan_theta_pos_lab>tanmax||tan_theta_ele_lab>tanmax))
	{
	  hCosThetaPosCoMPassing->Fill(cos_theta_pos_CoM);
	  reco_sqrts_straight_nores     = RecoSqrtsStraightTracks(pos_4mom_lab, ele_4mom_lab, 0);
	  reco_sqrts_straight_ECalSmear = RecoSqrtsStraightTracks(pos_4mom_lab, ele_4mom_lab, 1);
	}
      
      /*
      std::cout<<"Pos 1: ("<<x1_pos<<","<<y1_pos<<"), Ele 1: "<<x1_ele<<","<<y1_ele<<")"<<std::endl;
      std::cout<<"Pos 2: ("<<x2_pos<<","<<y2_pos<<"), Ele 2: "<<x2_ele<<","<<y2_ele<<")\n"<<std::endl;
      */
      /*      
      TVector3 pos_coord_nores1(x1_pos, y1_pos, z1);
      TVector3 pos_coord_nores2(x2_pos, y2_pos, z2);
      TVector3 ele_coord_nores1(x1_ele, y1_ele, z1);
      TVector3 ele_coord_nores2(x2_ele, y2_ele, z2);

      //simulate ECal resolution and add to e+/e- energy
      Double_t MaxECalRes_pos   = ECalResolution(Epos_true);
      Double_t MaxECalRes_ele   = ECalResolution(Eele_true);
      Double_t DeltaE_smear_pos = rand.Uniform(-1.*MaxECalRes_pos,MaxECalRes_pos);
      Double_t DeltaE_smear_ele = rand.Uniform(-1.*MaxECalRes_ele,MaxECalRes_ele);

      Double_t Epos_ECalSmear = Epos_true+DeltaE_smear_pos;
      Double_t Eele_ECalSmear = Eele_true+DeltaE_smear_ele;
      
      //s = (P1+P2)
      //no smearing
      Double_t reco_sqrts_nores = RecoSqrts(pos_coord_nores1, ele_coord_nores1, pos_coord_nores2, ele_coord_nores2, Epos_true, Eele_true);
      hDeltaSqrtsNoRes->Fill(reco_sqrts_nores-sqrts);

      //smearing just from ECal
      Double_t reco_sqrts_ECalSmear = RecoSqrts(pos_coord_nores1, ele_coord_nores1, pos_coord_nores2, ele_coord_nores2, Epos_ECalSmear, Eele_ECalSmear);
      hDeltaSqrtsECalSmear->Fill(reco_sqrts_ECalSmear-sqrts);

      //smearing from tracker
     */ 
    }

  TCanvas* cCosThetaPosCoMPassing = new TCanvas("cCosThetaPosCoMPassing","cCosThetaPosCoMPassing",900,700);
  hCosThetaPosCoMPassing->Draw();  

  TCanvas* cDeltaSqrtsNoRes = new TCanvas("cDeltaSqrtsNoRes","cDeltaSqrtsNoRes",900,700);
  hDeltaSqrtsNoRes->Draw();  

  TCanvas* cDeltaSqrtsECalSmear = new TCanvas("cDeltaSqrtsECalSmear","cDeltaSqrtsECalSmear",900,700);
  hDeltaSqrtsECalSmear->Draw();  
  
}
