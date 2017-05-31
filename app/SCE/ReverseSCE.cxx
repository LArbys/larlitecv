#include "ReverseSCE.h"

#include <string>

#include "TFile.h"

namespace larlitecv {

  ReverseSCE::ReverseSCE() {
    std::string sce_file_path = Form("%s/app/SCE/dat/reverse_sce_table.root",getenv("LARLITECV_BASEDIR"));
    rfile = new TFile(sce_file_path.c_str(),"OPEN");
    shift[0] = (TH3D*)rfile->Get("xorigin");
    shift[1] = (TH3D*)rfile->Get("yorigin");
    shift[2] = (TH3D*)rfile->Get("zorigin");
  }

  ReverseSCE::~ReverseSCE() {
    rfile->Close();
  }

  std::vector<double> ReverseSCE::getOffsets( const std::vector<double>& shiftedpos ) const {
    std::vector<double> originalpos = getOriginalPos( shiftedpos );
    for ( int i=0; i<3; i++ )
      originalpos[i] -= shiftedpos[i];
    return originalpos;
  }
  
  std::vector<float> ReverseSCE::getOffsets( const std::vector<float>& shiftedpos ) const {
    std::vector<float> originalpos = getOriginalPos( shiftedpos );
    for ( int i=0; i<3; i++ )
      originalpos[i] -= shiftedpos[i];;
    return originalpos;
  }
  
  std::vector<float> ReverseSCE::getOriginalPos( const std::vector<float>& shiftedpos ) const {
    std::vector<float> originalpos(3,0);
    int xbin = shift[0]->GetXaxis()->FindBin( shiftedpos[0] );
    int ybin = shift[0]->GetYaxis()->FindBin( shiftedpos[1] );
    int zbin = shift[0]->GetZaxis()->FindBin( shiftedpos[2] );

    if ( xbin==0 || xbin>=shift[0]->GetXaxis()->GetNbins()+1
	 || ybin==0 || ybin>=shift[0]->GetYaxis()->GetNbins()+1
	 || zbin==0 || zbin>=shift[0]->GetZaxis()->GetNbins()+1 ) {
      for (int i=0; i<3; i++)
	originalpos[i] = shiftedpos[i];
      return originalpos;
    }
    
    for (int i=0; i<3; i++)
      originalpos[i] = shiftedpos[i] + shift[i]->GetBinContent( xbin, ybin, zbin );

    return originalpos;    
  }
  
  std::vector<double> ReverseSCE::getOriginalPos( const std::vector<double>& shiftedpos ) const {
    std::vector<double> originalpos(3,0);
    int xbin = shift[0]->GetXaxis()->FindBin( shiftedpos[0] );
    int ybin = shift[0]->GetYaxis()->FindBin( shiftedpos[1] );
    int zbin = shift[0]->GetZaxis()->FindBin( shiftedpos[2] );

    if ( xbin==0 || xbin>=shift[0]->GetXaxis()->GetNbins()+1
	 || ybin==0 || ybin>=shift[0]->GetYaxis()->GetNbins()+1
	 || zbin==0 || zbin>=shift[0]->GetZaxis()->GetNbins()+1 ) {
      for (int i=0; i<3; i++)
	originalpos[i] = shiftedpos[i];
      return originalpos;
    }
    
    for (int i=0; i<3; i++)
      originalpos[i] = shiftedpos[i] + shift[i]->GetBinContent( xbin, ybin, zbin );

    return originalpos;
  }

}
