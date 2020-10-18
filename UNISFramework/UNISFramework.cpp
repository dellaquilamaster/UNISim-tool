#include <UNISFramework.h>

//____________________________________________________
UNISFramework::UNISFramework() : UNISFramework("config/UNISim.conf")
{}

//____________________________________________________
UNISFramework::UNISFramework(const char * file_config) :
fConfigurationFile(file_config),
fBeamCenter(0,0,0),
fBeamAngularSpread(0),
fBeamPositionSpread(0),
fPhysicsModelName(""),
fPhysicsConfigFileName(""),
fOutputFolder("./output/"),
fOutputFileName(),
fGraphics(false),
fAdvancedGraphics(false),
fTheTree(0),
fTheFile(0),
fUMAToMeV(931.4936148),
fcmToum(10000),
fTheEventGenerator(0),
fApp(new TRint("The Unified-Simulation-tool", new int(), 0, 0, 0, kTRUE))
{
  gNucData = new nuclear_masses("./LISETools/input/masses.conf");
  gLISEELossModule = new EnergyLossModule("./LISETools/");
}

//____________________________________________________
UNISFramework::~UNISFramework()
{
  fApp->Terminate(0);
}

//____________________________________________________
void UNISFramework::SetIterations(int value)
{
  fNumEvents=value;
}

//____________________________________________________
int UNISFramework::ConfigureFramework()
{
  std::ifstream FileIn(fConfigurationFile.c_str());
  if(!FileIn.is_open()) {
    return -1;
  }
  
  int NRead=0;

  while (!FileIn.eof())
  {
    std::string LineRead;
    std::getline(FileIn, LineRead);

    std::string LineReadCommentLess (LineRead.substr(0, LineRead.find("*")));

    if(LineReadCommentLess.empty()) continue;
    if(LineReadCommentLess.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineReadCommentLess);
    
    std::string Command;
    
    LineStream >> Command;
        
    if(Command.compare("set")==0) {
      NRead+=ProcessSetCommand(LineReadCommentLess.c_str());
    } else if (Command.compare("define")==0) {
      NRead+=ProcessDefineCommand(LineReadCommentLess.c_str());
    } else if (Command.compare("add")==0) {
      NRead+=ProcessAddCommand(LineReadCommentLess.c_str());
    }
  }
  
  return NRead;  
}

//____________________________________________________
int UNISFramework::ProcessSetCommand(const char * line)
{
  std::string InputLine(line);
  std::istringstream LineStream(InputLine);
  
  std::string Command;
  std::string WhatToSet;
  
  LineStream >> Command >> WhatToSet;
  
  std::string ValueToSet;
  
  LineStream>>ValueToSet;
  
  if(WhatToSet.compare("VERBOSE_MODE")==0) {
    fVerbose=ValueToSet.compare("true")==0 ? true : false;
  } else if(WhatToSet.compare("GRAPHICAL_MODE")==0) {
    if(ValueToSet.compare("light")==0) {
      fGraphics=true; 
    } else if(ValueToSet.compare("advanced")==0) {
      fGraphics=true;
      fAdvancedGraphics=true;
    }
    fGraphics=ValueToSet.compare("false")==0 ? false : true;
  } else if(WhatToSet.compare("OUTPUT_DIRECTORY")==0) {
    fOutputFolder.assign(ValueToSet.substr(ValueToSet.find("\"")+1,ValueToSet.find_last_of("\"")-(ValueToSet.find("\"")+1)));
    if(fOutputFolder.find_last_of('/')!=fOutputFolder.length()-1) {
      fOutputFolder.append("/");
    }
  } else if(WhatToSet.compare("RANDOM_SEED")==0) {
    gRandomSeed=std::stof(ValueToSet);
    gRandom->SetSeed(gRandomSeed);
  } else if(WhatToSet.compare("BEAM")==0) {
    do {
      if(ValueToSet.find("-Z")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-Z=")+3)); 
        fTheBeam.fZ=std::stoi(ValueToSet); 
      } else if(ValueToSet.find("-A")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-A=")+3));
        fTheBeam.fA=std::stoi(ValueToSet); 
      }
    } while (LineStream>>ValueToSet);
  } else if(WhatToSet.compare("TARGET")==0) {
    do {
      if(ValueToSet.find("-Z")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-Z=")+3)); 
        fTheTarget.fZ=std::stoi(ValueToSet); 
      } else if(ValueToSet.find("-A")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-A=")+3));
        fTheTarget.fA=std::stoi(ValueToSet); 
      }
    } while (LineStream>>ValueToSet);
  } else if(WhatToSet.compare("BEAM_ENERGY")==0) {
    fBeamEnergy=std::stof(ValueToSet);
  } else if(WhatToSet.compare("BEAM_POSITION")==0) {
    do {
      if(ValueToSet.find("-X")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-X=")+3)); 
        fBeamCenter.SetX(std::stof(ValueToSet)); 
      } else if(ValueToSet.find("-Y")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-Y=")+3));
        fBeamCenter.SetY(std::stof(ValueToSet)); 
      } else if(ValueToSet.find("-Z")!=std::string::npos) {
        ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-Z=")+3));
        fBeamCenter.SetZ(std::stof(ValueToSet)); 
      }
    } while (LineStream>>ValueToSet);
  } else if(WhatToSet.compare("BEAM_ANGULAR_SPREAD")==0) {
    fBeamAngularSpread=std::stof(ValueToSet)*TMath::DegToRad();
  } else if(WhatToSet.compare("BEAM_POSITION_SPREAD")==0) {
    fBeamPositionSpread=std::stof(ValueToSet);
  } else if(WhatToSet.compare("TARGET_MATERIAL")==0) {
    fTargetMaterial.assign(ValueToSet);
  } else if(WhatToSet.compare("TARGET_THICKNESS")==0) {
    fTargetThickness=std::stof(ValueToSet);
  } else if(WhatToSet.compare("PHYSICS_MODEL")==0) {
      if( (fPhysicsModelName.empty() && ValueToSet.compare("SequentialDecay")==0) || fPhysicsModelName.compare("SequentialDecay")==0) {
        fPhysicsModelName.assign("SequentialDecay");
        fTheEventGenerator = new UNISSequentialDecay();
        LineStream>>ValueToSet;
        if(fPhysicsConfigFileName.empty()) fPhysicsConfigFileName.assign(ValueToSet.substr(ValueToSet.find("\"")+1,ValueToSet.find_last_of("\"")-(ValueToSet.find("\"")+1)));
        if(fTheEventGenerator->LoadConfiguration(fPhysicsConfigFileName.c_str())<=0) {
          printf("Error: error while building SequentialDecay event generator from file %s\nAborting!", fPhysicsConfigFileName.c_str());
          exit(1);
        }
      } else if( (fPhysicsModelName.empty() && ValueToSet.compare("RutherfordScattering")==0) || fPhysicsModelName.compare("RutherfordScattering")==0) {
        fPhysicsModelName.assign("RutherfordScattering");
        fTheEventGenerator = new UNISRutherfordScattering();
        LineStream>>ValueToSet;
        if(fPhysicsConfigFileName.empty()) fPhysicsConfigFileName.assign(ValueToSet.substr(ValueToSet.find("\"")+1,ValueToSet.find_last_of("\"")-(ValueToSet.find("\"")+1)));
        if(fTheEventGenerator->LoadConfiguration(fPhysicsConfigFileName.c_str())<=0) {
          printf("Error: error while building RutherfordScattering event generator from file %s\nAborting!", fPhysicsConfigFileName.c_str());
          exit(1);
        }
      } else if( (fPhysicsModelName.empty() && ValueToSet.compare("SequentialDecayTwoBody")==0) || fPhysicsModelName.compare("SequentialDecayTwoBody")==0) {
        fPhysicsModelName.assign("SequentialDecayTwoBody");
        fTheEventGenerator = new UNISSequentialDecayTwoBody();
        LineStream>>ValueToSet;
        if(fPhysicsConfigFileName.empty()) fPhysicsConfigFileName.assign(ValueToSet.substr(ValueToSet.find("\"")+1,ValueToSet.find_last_of("\"")-(ValueToSet.find("\"")+1)));
        if(fTheEventGenerator->LoadConfiguration(fPhysicsConfigFileName.c_str())<=0) {
          printf("Error: error while building SequentialDecayTwoBody event generator from file %s\nAborting!", fPhysicsConfigFileName.c_str());
          exit(1);
        }
      } else return 0;
  } else {
    return 0; 
  }
  
  return 1;
}

//____________________________________________________
int UNISFramework::ProcessDefineCommand(const char * line)
{
  std::string InputLine(line);
  std::istringstream LineStream(InputLine);  
  
  std::string Command;
  std::string WhatToSet;
  
  LineStream>>Command>>WhatToSet;
  
  if(WhatToSet.compare("SETUP")==0) {
    std::string SetupName;
    LineStream>>SetupName;
    SetupName.assign(SetupName.substr(SetupName.find("\"")+1,SetupName.find_last_of("\"")-(SetupName.find("\"")+1)));
    
    UNISDetectionSetup * NewSetup = new UNISDetectionSetup(SetupName.c_str());
    
    fExpSetup=NewSetup;
    
  } else {
    return 0; 
  }
  
  return 1;
}

//____________________________________________________
int UNISFramework::ProcessAddCommand(const char * line)
{
  std::string InputLine(line);
  std::istringstream LineStream(InputLine);
  
  std::string Command;
  std::string WhatToSet;
  
  LineStream>>Command>>WhatToSet;
  
  if(WhatToSet.compare("DETECTOR")==0) {
    //Adding a new detector
    std::string DetectorType;
    LineStream>>DetectorType;
    std::string ValueToSet;
        
    if(DetectorType.compare("DSSSD")==0) {
      double distance=0;
      double theta_pos=0;
      double phi_pos=0;
      int strip_number=0;
      double strip_width=0;
      double strip_inter=0;
      double frame_width=0;
      double dead_layer=0;
      
      while (LineStream>>ValueToSet) {
        if(ValueToSet.find("-distance")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-distance=")+10)); 
          distance=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-theta")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-theta=")+7)); 
          theta_pos=std::stof(ValueToSet)*TMath::DegToRad();
        } else if(ValueToSet.find("-phi")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-phi=")+5)); 
          phi_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-strips")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-strips=")+8)); 
          strip_number=std::stoi(ValueToSet); 
        } else if(ValueToSet.find("-strip_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-strip_width=")+13)); 
          strip_width=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-inter_strip")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-inter_strip=")+13)); 
          strip_inter=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-frame_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-frame_width=")+13)); 
          frame_width=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-dead_layer")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-dead_layer=")+12)); 
          dead_layer=std::stof(ValueToSet); 
        }
      }
      
      UNISStripDetector * NewDetector = new UNISStripDetector(distance,theta_pos,phi_pos,strip_number,strip_width,strip_inter,frame_width,dead_layer,"");
      fExpSetup->RegisterUnit(NewDetector);
      
    } else if(DetectorType.compare("DSSSD_ROT")==0) {
      double X0=0;
      double Y0=0;
      double Z0=0;
      double tilt_X=0;
      double tilt_Y=0;
      int strip_number=0;
      double strip_width=0;
      double strip_inter=0;
      double frame_width=0;
      double dead_layer=0;
      
      while (LineStream>>ValueToSet) {
        if(ValueToSet.find("-X0")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-X0=")+4)); 
          X0=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-Y0")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-Y0=")+4)); 
          Y0=std::stof(ValueToSet);
        } else if(ValueToSet.find("-Z0")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-Z0=")+4)); 
          Z0=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-tilt_X")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-tilt_X=")+8)); 
          tilt_X=std::stof(ValueToSet)*TMath::DegToRad();
        } else if(ValueToSet.find("-tilt_Y")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-tilt_Y=")+8)); 
          tilt_Y=std::stof(ValueToSet)*TMath::DegToRad();
        } else if(ValueToSet.find("-strips")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-strips=")+8)); 
          strip_number=std::stoi(ValueToSet); 
        } else if(ValueToSet.find("-strip_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-strip_width=")+13)); 
          strip_width=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-inter_strip")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-inter_strip=")+13)); 
          strip_inter=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-frame_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-frame_width=")+13)); 
          frame_width=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-dead_layer")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-dead_layer=")+12)); 
          dead_layer=std::stof(ValueToSet); 
        }
      }
      
      UNISStripDetector * NewDetector = new UNISStripDetector(X0,Y0,Z0,tilt_X,tilt_Y,strip_number,strip_width,strip_inter,frame_width,dead_layer,"");
      fExpSetup->RegisterUnit(NewDetector);
      
    } else if(DetectorType.compare("LAMP_WEDGE")==0) {
      double distance=0;
      double phi_pos=0;
      double tilt=0;
      double bottom_frame_distance=0;
      int strip_number=0;
      double strip_width=0;
      double strip_inter=0;
      while (LineStream>>ValueToSet) {
        if(ValueToSet.find("-distance")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-distance=")+10)); 
          distance=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-phi_pos")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-phi_pos=")+9)); 
          phi_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-tilt")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-tilt=")+6)); 
          tilt=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-frame_distance")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-frame_distance=")+16)); 
          bottom_frame_distance=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-strips")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-strips=")+8)); 
          strip_number=std::stoi(ValueToSet); 
        } else if(ValueToSet.find("-strip_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-strip_width=")+13)); 
          strip_width=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-inter_strip")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-inter_strip=")+13)); 
          strip_inter=std::stof(ValueToSet); 
        }
      }
      
      UNISLampWedgeDetector * NewDetector = new UNISLampWedgeDetector(distance,phi_pos,tilt,bottom_frame_distance,strip_number,strip_width,strip_inter);
      fExpSetup->RegisterUnit(NewDetector);
    } else if(DetectorType.compare("LAMP_WEDGE_MMM")==0) {
      double distance=0;
      double phi_pos=0;
      double tilt=0;
      double bottom_frame_distance=0;
      double strip_inter=0;
      while (LineStream>>ValueToSet) {
        if(ValueToSet.find("-distance")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-distance=")+10)); 
          distance=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-phi_pos")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-phi_pos=")+9)); 
          phi_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-tilt")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-tilt=")+6)); 
          tilt=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-frame_distance")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-frame_distance=")+16)); 
          bottom_frame_distance=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-inter_strip")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-inter_strip=")+13)); 
          strip_inter=std::stof(ValueToSet); 
        }
      }
      
      UNISLampWedgeMMMDetector * NewDetector = new UNISLampWedgeMMMDetector(distance,phi_pos,tilt,bottom_frame_distance,strip_inter);
      fExpSetup->RegisterUnit(NewDetector);
    } else if(DetectorType.compare("FAZIA_BLOCK")==0) {
      double displacement=0;
      double theta_pos=0;
      double phi_pos=0;
      double pad_width=0;
      double frame_width=0;
      while (LineStream>>ValueToSet) {
        if(ValueToSet.find("-displacement")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-displacement=")+14));
          displacement=std::stoi(ValueToSet); 
        } else if(ValueToSet.find("-phi_pos")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-phi_pos=")+9)); 
          phi_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-theta_pos")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-theta_pos=")+11)); 
          theta_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-pad_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-pad_width=")+11)); 
          pad_width=std::stof(ValueToSet); 
        } else if(ValueToSet.find("-frame_width")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-frame_width=")+13)); 
          frame_width=std::stof(ValueToSet); 
        }
      }
      UNISFaziaBlock * NewDetector = new UNISFaziaBlock(theta_pos,phi_pos,displacement,pad_width,frame_width);
      fExpSetup->RegisterUnit(NewDetector);
    } else if(DetectorType.compare("PHOTO_DIODE")==0) {
      double distance=0;
      double theta_pos=0;
      double phi_pos=0;
      while (LineStream>>ValueToSet) {
        if(ValueToSet.find("-distance")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-distance=")+10));
          distance=std::stoi(ValueToSet); 
        } else if(ValueToSet.find("-phi_pos")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-phi_pos=")+9)); 
          phi_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } else if(ValueToSet.find("-theta_pos")!=std::string::npos) {
          ValueToSet.assign(ValueToSet.substr(ValueToSet.find("-theta_pos=")+11)); 
          theta_pos=std::stof(ValueToSet)*TMath::DegToRad(); 
        } 
      }
      UNISSiliconPhotoDiode * NewDetector = new UNISSiliconPhotoDiode(theta_pos,phi_pos,distance);
      fExpSetup->RegisterUnit(NewDetector);
    }
    
  } else {
    return 0; 
  }
  
  return 1;
}

//____________________________________________________
void UNISFramework::PrintConfiguration() const
{
  printf("\nThe UNISim-tool Framework\n");
  
  //
  //Printing date and time
  time_t now;
  time(&now);
  printf("Date and time: %s\n", ctime(&now));
  //
  
  //
  //Printing configuration quantities
  printf("----------------------\n");
  printf("Verbose Mode: %s\n", fVerbose ? "on" : "off");
  printf("Graphical Mode: %s\n", fGraphics&&!fAdvancedGraphics ? "on (light)" : (fGraphics&&fAdvancedGraphics ? "on (advanced)" : "off"));
  printf("----------------------\n");
  printf("Beam: Z=%d A=%d Ekin=%f\n", fTheBeam.fZ, fTheBeam.fA, fBeamEnergy);
  printf("Target: Z=%d A=%d\n", fTheTarget.fZ, fTheTarget.fA);
  printf("Physics Model: %s\n", fPhysicsModelName.c_str());
  printf("Physics Configuration: %s\n", fPhysicsConfigFileName.c_str());
  printf("****\n");
  printf("----------------------\n");
  printf("****\n");
  printf("----------------------\n");
  printf("Detection Setup: \"%s\" (%d detectors)\n", fExpSetup->GetName(), fExpSetup->Size());
  printf("----------------------\n");
  //
  
  printf("\n");
}

//____________________________________________________
void UNISFramework::InitTree()
{
  //
  //Creating output file
  if(fOutputFileName.empty()) fOutputFileName.assign(Form("UNIS_%devents.root",fNumEvents));
  fTheFile = new TFile(Form("%s%s",fOutputFolder.c_str(),fOutputFileName.c_str()),"RECREATE");
  //
  
  //
  //Allocating output tree
  fTheTree = new TTree("t", "UNISim-tool Framework Data");
  //
  
  //
  //Allocating data structure
  fevt = new UNISRootEvent();
  //
  
  //
  //Creating Thee branch  
  fTheTree->Branch(Form("%s.", fExpSetup->GetName()), "UNISRootEvent", fevt, 3200, 2);
  //
  
  //
  //Setting Autosave
  fTheTree->SetAutoSave(500000);
  //
}

//____________________________________________________
void UNISFramework::EndProcess()
{
  //
  //Saving tree to file
  fTheTree->AutoSave();
  //
  
  //
  //Closing Root File
  fTheFile->Close();
  //
  
  //
  printf("\n\nSuccessfully closed file %s\n\n", Form("%s%s",fOutputFolder.c_str(),fOutputFileName.c_str()));
  //
}

//____________________________________________________
int UNISFramework::ReadInput(int argc, char ** argv)
{
  int read=0;
  for(int i=1; i<argc; i++) {
    if(strcmp(argv[i],"-events")==0) {
      SetIterations(atoi(argv[++i]));
      read++;
    }
    if(strcmp(argv[i],"-o")==0) {
      fOutputFileName.assign(argv[++i]);
      read++;
    }
    if(strcmp(argv[i],"-physics")==0) {
      fPhysicsModelName.assign(argv[++i]);
      read++;
    }
    if(strcmp(argv[i],"-reaction")==0) {
      fPhysicsConfigFileName.assign(argv[++i]);
      read++;
    }
  }
  
  return read;
}

//____________________________________________________
void UNISFramework::ProcessIterations()
{
  //This is the core method.
  //Here we process, event-by-event, the framework iterations.
  //The tree is filled for each event.
  
  //
  //Initialization of the target
  InitializeTarget();
  //
  
  if(fGraphics) {
    fExpSetup->Draw3D(fAdvancedGraphics ? "ogl" : "");
  }
  
  //
  fStartTime=clock();
  //
  
  //
  //Main Iteration
  for(int ientry=0; ientry<fNumEvents; ientry++)
  {
    //
    //
    //
    
    //
    //Printing progress bar
    if(ientry%100==0) {
      PrintPercentage(ientry, fNumEvents);
    }
    //
    
    //
    //Verbose mode
    if(fVerbose) {
      printf("*** Beginning of iteration %d\n", ientry);
      printf("**************************************************\n\n");
    }
    //
    
    //
    //Setting the projectile for the iteration
    GenerateBeam(); //Initialization of the beam particle producing the reaction (event-by-event)
    //
    
    //
    //Generate Event
    std::vector<UNISIon> TheEvent = fTheEventGenerator->GetEvent();
    //
    
    //
    //Writing the event on the UNISRootEvent data structure
    RegisterEvent(TheEvent);
    //
    
    //
    //Detect emitted particles
    DetectEvent();
    //
    
    //
    //Display Tracks if needed
    if(fGraphics && !fAdvancedGraphics) {  
      DisplayTracks(); 
    }
    //
    
    //
    //Fill the output tree
    FillEvent();
    //
  }
  //End of the Main Iteration
  
  //
  //Running TEveManager standalone
  if(fGraphics) {
    fApp->Run(kTRUE);
    fApp->Terminate(0);
  }
  //
}

//____________________________________________________
void UNISFramework::GenerateBeam()
{
  //
  //Creation of the beam
  fTheBeam.fMass = gNucData->get_mass_Z_A(fTheBeam.fZ,fTheBeam.fA);
  //Interaction with the target
  fBeamEnergyMidTarget = fTargetThickness>0 ? fBeamEnergy - gLISEELossModule->GetEnergyLoss(fTheBeam.fZ,fTheBeam.fA,fBeamEnergy,fTargetMaterial.c_str(),fTargetThickness) : fBeamEnergy;
  double BeamMomentum = sqrt(pow(fBeamEnergyMidTarget+fTheBeam.fMass,2)-pow(fTheBeam.fMass,2));
  fTheBeam.fMomentum=TLorentzVector(0,0,BeamMomentum,fBeamEnergyMidTarget+fTheBeam.fMass);
  fTheBeam.fMomentum.RotateY(gRandom->Gaus(0.,fBeamAngularSpread/2.355)); //beam deviation (polar angle)
  fTheBeam.fMomentum.RotateZ(gRandom->Uniform(0.,2*TMath::Pi())); //randomization of the beam direction around the Z axis (cylindrical symmetry)
  fBeamPosition.SetXYZ(fBeamCenter.X()+gRandom->Gaus(0,fBeamPositionSpread/2.355),fBeamCenter.Y()+gRandom->Gaus(0,fBeamPositionSpread/2.355),fBeamCenter.Z());
  //
  //Initialization of the beam for the simulation
  fTheEventGenerator->SetBeam(fTheBeam);
  //
  
}

//____________________________________________________
void UNISFramework::InitializeTarget()
{
  //
  //Creation of the target
  fTheTarget.fMass = gNucData->get_mass_Z_A(fTheTarget.fZ,fTheTarget.fA);
  fTheTarget.fMomentum=TLorentzVector(0,0,0,fTheTarget.fMass);
  //
  //Initialization of the target
  fTheEventGenerator->SetTarget(fTheTarget);
  //
}

//____________________________________________________
void UNISFramework::DetectEvent()
{
  //
  //Loop on the registered event to process the detection
  for(int i=0; i<fevt->fmulti; i++) {
    int det_index = fExpSetup->GetDetectorIndex(fevt->fThetaOrigin[i],fevt->fPhiOrigin[i],fBeamPosition.X(),fBeamPosition.Y(),fBeamPosition.Z());
    if(det_index>=0) {
      UNISDetectionUnit * TheDetector = fExpSetup->GetDetector(det_index);
      int num_pixel = TheDetector->GetPixel(fevt->fThetaOrigin[i],fevt->fPhiOrigin[i],fBeamPosition.X(),fBeamPosition.Y(),fBeamPosition.Z());
      if(num_pixel>=0) {
        //The particle is detected
        fevt->fmulti_detected++;
        fevt->fIsDetected[i]=true;
        TVector3 impact_point = TheDetector->GetImpactPointLab(fevt->fThetaOrigin[i],fevt->fPhiOrigin[i],fBeamPosition.X(),fBeamPosition.Y(),fBeamPosition.Z());
        fevt->fnumdet[i]=det_index;
        fevt->fnumpixel[i]=num_pixel;
        fevt->fXDetHit[i]=impact_point.X();
        fevt->fYDetHit[i]=impact_point.Y();
        fevt->fZDetHit[i]=impact_point.Z();
      }
    }
  }
}

//____________________________________________________
void UNISFramework::RegisterEvent(std::vector<UNISIon> & AnEvent)
{
  fevt->fmulti=0;
  fevt->fmulti_detected=0;
  
  //
  //Loop on the event
  for(int i=0; i<AnEvent.size(); i++) {       
    fevt->fIsDetected[i]=false;
    fevt->fnumdet[i]=-1;
    fevt->fnumpixel[i]=-1;
    fevt->fKinEnergy[i]=AnEvent[i].fMomentum.E()-AnEvent[i].fMomentum.M();
    fevt->fXDetHit[i]=-9999;
    fevt->fYDetHit[i]=-9999;
    fevt->fZDetHit[i]=-9999;
    fevt->fKinEnergyOrigin[i]=AnEvent[i].fMomentum.E()-AnEvent[i].fMomentum.M();
    fevt->fThetaOrigin[i]=AnEvent[i].fMomentum.Theta();
    fevt->fPhiOrigin[i]=AnEvent[i].fMomentum.Phi();
    fevt->fZ[i]=AnEvent[i].fZ;
    fevt->fA[i]=AnEvent[i].fA;
    fevt->fKinEnergyAfterTarget[i]=(fTargetThickness>0 ? (fevt->fKinEnergyOrigin[i] - (cos(fevt->fThetaOrigin[i])!=0 ? gLISEELossModule->GetEnergyLoss(fevt->fZ[i],fevt->fA[i],fevt->fKinEnergyOrigin[i],fTargetMaterial.c_str(),fTargetThickness/2./cos(fevt->fThetaOrigin[i])) : fevt->fKinEnergyOrigin[i])) : fevt->fKinEnergyOrigin[i]);
    fevt->fKinEnergyOriginCms[i]=-9999;
    fevt->fThetaOriginCms[i]=-9999;
    fevt->fmulti++;
  }
  //
}

//____________________________________________________
void UNISFramework::FillEvent()
{
  //
  //Fill output tree
  fTheTree->Fill();
  //
}

//____________________________________________________
void UNISFramework::DisplayTracks()
{
  //
  //Display the beam track
  TGeoTrack * BeamTrack = new TGeoTrack();
  BeamTrack->AddPoint(-fBeamCenter.Y(),fBeamCenter.X(),-20,0);
  BeamTrack->AddPoint(-fBeamCenter.Y(),fBeamCenter.X(),fBeamCenter.Z(),0);
  BeamTrack->SetLineColor(kBlack);
  gGeoManager->AddTrack(BeamTrack);
  //
  
  //
  //Looping on event multiplicity to draw tracks
  for(int i=0; i<fevt->fmulti; i++) {
    TGeoTrack * ParticleTrack = new TGeoTrack();
    if(fevt->fXDetHit[i]!=-9999 && fevt->fYDetHit[i]!=-9999 && fevt->fZDetHit[i]!=-9999) { //The particle is detected
      ParticleTrack->AddPoint(-fBeamCenter.Y(),fBeamCenter.X(),fBeamCenter.Z(),0);
      ParticleTrack->AddPoint(-fevt->fYDetHit[i],fevt->fXDetHit[i],fevt->fZDetHit[i],0);
    } else { //The particle escaped from the region
      ParticleTrack->AddPoint(-fBeamCenter.Y(),fBeamCenter.X(),fBeamCenter.Z(),0);
      ParticleTrack->AddPoint(-20*sin(fevt->fThetaOrigin[i])*sin(fevt->fPhiOrigin[i]),20*sin(fevt->fThetaOrigin[i])*cos(fevt->fPhiOrigin[i]),20*cos(fevt->fThetaOrigin[i]),0);
    }
    fevt->fZ[i]<=10 ? ParticleTrack->SetLineColor(fevt->fZ[i]+fevt->fA[i]) : ParticleTrack->SetLineColor(kRed);
    gGeoManager->AddTrack(ParticleTrack);
  }
  //
  
  gGeoManager->DrawTracks();
}
  
//____________________________________________________
void UNISFramework::PrintPercentage(int jentry, int nentries) const
{
  double time_elapsed = (double)(clock() - fStartTime)/CLOCKS_PER_SEC;
  std::cout << "  Percentage = " << std::fixed << std::setprecision(1) << std::setw(5) << (100*double(jentry)/nentries) << " %";
  std::cout << "   [";
  int printindex=0;
  for(; printindex<int(100*double(jentry)/nentries); printindex+=5) std::cout << "=";
  for(; printindex<100; printindex+=5) std::cout << " ";
  std::cout << "]   " << "elapsed time " << std::setprecision(1) <<
  (time_elapsed<60 ? time_elapsed : (time_elapsed<3600 ? time_elapsed/60 : time_elapsed/3600)) <<
  (time_elapsed<60 ? " s; " : (time_elapsed<3600 ? " m; " : " h; "));
  if(jentry>0) {
    double time_remaining = (time_elapsed/jentry)*(nentries-jentry);
    std::cout << " estimated remaining time " << std::setprecision(1) <<
    (time_remaining<60 ? time_remaining : (time_remaining<3600 ? time_remaining/60 : time_remaining/3600)) <<
    (time_remaining<60 ? " s      " : (time_remaining<3600 ? " m      " : " h      "));
  }
  std::cout << "\r";
  std::cout.flush();
}
