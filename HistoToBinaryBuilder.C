void HistoToBinaryBuilder()
{
  //file containing 2D histo of bremsstrahlung energy/angle
  TString sBremFilePath = "/Users/elizabeth/Documents/PADME/ThesisRootFiles/Quick_MC_Run660_AnalysisOutputFile.root";
  TFile* fBremFile = TFile::Open(sBremFilePath);
  if(!fBremFile || fBremFile->IsZombie())
    {
      std::cerr << "[BremArrayBuilder] No file " << sBremFilePath << std::endl;
      return;
    }

  //path of 2D histo containing bremsstrahlung energy/angle distribution
  TString sHistoPath = "MCTruth/Bremstrahlung EvsAn";
  TH2* hBremEvsAn = (TH2*)fBremFile->Get(sHistoPath);
  if (!hBremEvsAn) {
    std::cout<<"[BremArrayBuilder] No histogram "<<sHistoPath<<" in "<<sBremFilePath<<std::endl;
    fBremFile->Close(); delete fBremFile;
    return;
  }

  TCanvas* cBremEvsAn = new TCanvas("cBremVsAn","cBremVsAn",900,700);
  hBremEvsAn->Draw("colz");
  cBremEvsAn->SetLogz();

  const int nAngleBins = hBremEvsAn->GetNbinsX();
  const int nEnergyBins = hBremEvsAn->GetNbinsY();

  // Build CDF (flattened)
  std::vector<double> vCDF;
  vCDF.reserve(nAngleBins * nEnergyBins);

  double cumulative = 0.;
  for(int iEnergy = 1; iEnergy<=nEnergyBins; iEnergy++)
    {
      for(int iAngle = 1; iAngle<=nAngleBins; iAngle++)
	{
	  cumulative += hBremEvsAn->GetBinContent(iAngle, iEnergy);

	  //flattened index = iEnergy*nAngleBins + iAngle;
	  vCDF.push_back(cumulative);
	}
    }

  // Extract bin edges
  std::vector<double> AngleBinEdges(nAngleBins + 1);
  std::vector<double> EnergyBinEdges(nEnergyBins + 1);

  for (int iAngle = 1; iAngle <= nAngleBins + 1; ++iAngle)
    AngleBinEdges[iAngle - 1] = hBremEvsAn->GetXaxis()->GetBinLowEdge(iAngle);

  for (int iEnergy = 1; iEnergy <= nEnergyBins + 1; ++iEnergy)
    EnergyBinEdges[iEnergy - 1] = hBremEvsAn->GetYaxis()->GetBinLowEdge(iEnergy);

  // Write binary file
  std::ofstream out("BremSampler.bin", std::ios::binary);
  if (!out.is_open())
  {
    std::cerr << "[BremArrayBuilder] Failed to open output file\n";
    return;
  }

  // Write dimensions
  out.write(reinterpret_cast<const char*>(&nAngleBins), sizeof(int));
  out.write(reinterpret_cast<const char*>(&nEnergyBins), sizeof(int));

  // Write x edges
  out.write(reinterpret_cast<const char*>(AngleBinEdges.data()),
            AngleBinEdges.size() * sizeof(double));

  // Write y edges
  out.write(reinterpret_cast<const char*>(EnergyBinEdges.data()),
            EnergyBinEdges.size() * sizeof(double));

  // Write CDF
  out.write(reinterpret_cast<const char*>(vCDF.data()),
            vCDF.size() * sizeof(double));

  out.close();

  std::cout << "[BremArrayBuilder] Written BremSampler.bin\n";
}
