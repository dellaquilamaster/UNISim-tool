#include "UNISTarget.h"
#include "UNISTargetStack.h"
#include "TFile.h"

int main() {
    gLISEELossModule = new EnergyLossModule("../../LISETools/input/");
    //gLISEELossModule->SetEnergyLossPrecision(0.001);
    gNucData = new nuclear_masses("../../LISETools/input/masses.conf");

    //Target target;
    //target.material = "C";
    //target.IsActive = true;
    //target.thickness_um = 0.5;
    //target.thickness_ug = 500;

    //printf("target.material: %s\n", target.material.c_str());
    //printf("target.IsActive: %d\n", target.IsActive);
    //printf("target.thickness_um: %f\n", target.thickness_um);
    //printf("target.thickness_ug: %f\n", target.thickness_ug);

    //UNISTarget t;
    //t.SetName("Carbon_target");  
    //t.SetMaterial("12C");
    //t.SetThickness(1.5, 1000.0);
    //t.SetIsActive(true);
    //t.Print();

    //if(t.TargetIsActive())
    //    printf("Aktivna\n");
    //else
    //    printf("Nije Aktivna\n");

    //t.GetTargetMaterial();
    //t.GetTargetIsActive();
    //t.GetTargetThickness_um();
    //t.GetTargetThickness_ug();
  
    gRandomSeed = 10;
    UNISTargetStack targetStack_1;
    targetStack_1.SetName("SkupMeta_1");

    UNISTarget t1("prva_meta", "Al", false, 3.330866, 900);
    UNISTarget t2("druga_meta", "B", true, 4.254957, 1000);

    targetStack_1.AddTarget(&t1);
    targetStack_1.AddTarget(&t2);

//______________________Test Point Of Interaction___________________________________________________

    //This stack is perfect to test interaction point production
    UNISTargetStack testStack;
    Interaction_Point p0 = testStack.GetInteractionPoint();
    printf("p0.target_number = %d\np0.depth = %f\n", p0.target_number, p0.depth);

    UNISTargetStack targetStack_2;
    targetStack_2.SetName("SkupMeta_2");
//
    UNISTarget ts_1("first_meta", "B", true, 1.0637393, 250);
    UNISTarget ts_2("second_meta", "Al", false, 0.14803849, 40);
    UNISTarget ts_3("third_meta", "B", true, 0.4254957, 100);
    UNISTarget ts_4("fourth_meta", "Al", false, 0.14803849, 40);
    UNISTarget ts_5("fifth_meta", "B", true, 2.1274785, 500);
    UNISTarget ts_6("sixth_meta", "Al", false, 0.14803849, 40);
    UNISTarget ts_7("seventh_meta", "B", true, 1.0637393, 250);
    UNISTarget ts_8("eigth_meta", "Al", false, 0.14803849, 40);
//
    targetStack_2.AddTarget(&ts_1);
    targetStack_2.AddTarget(&ts_2);
    targetStack_2.AddTarget(&ts_3); 
    targetStack_2.AddTarget(&ts_4); 
    targetStack_2.AddTarget(&ts_5);
    targetStack_2.AddTarget(&ts_6);
    targetStack_2.AddTarget(&ts_7);
    targetStack_2.AddTarget(&ts_8);

    //targetStack_2.RemoveTarget(&ts_2);
    //targetStack_2.ListStack();
    //targetStack_2.RemoveTarget(&ts_4);
    //targetStack_2.ListStack();
    //targetStack_2.RemoveTarget(&ts_6);
    //targetStack_2.ListStack();
    //targetStack_2.RemoveTarget(&ts_8);
    //targetStack_2.ListStack();

    TFile* outFile = new TFile("./test.root", "RECREATE");
    double totalThickenss_um = 0;
    for(int i = 0; i<targetStack_2.GetNum_ofTargets_inStack(); ++i)
    {
        totalThickenss_um += targetStack_2.GetTarget_fromStack(i).GetTargetThickness_um();
    }
    TH1D* count_per_target = new TH1D("count_per_target", "count_per_target; target; count", 8, -0.5, 7.5);
    TH1D* count_per_distance = new TH1D("count_per_distance", "count_per_distance; dist_in_um; count", 100, 0, totalThickenss_um);

    Interaction_Point p;
    for(int i = 0; i< 110000; i++)
    {
        p = targetStack_2.GetInteractionPoint();
        count_per_target->Fill(p.target_number);
        count_per_distance->Fill(p.depth);
    }

    outFile->cd();
    count_per_target->Write();
    count_per_distance->Write();
    outFile->Close();
//___________________________________________________________________________


    //Add, List, Remove, GetNormalVector, SetOrientation ✓
    //targetStack_1.ListStack();    
    //targetStack_2.ListStack();    

    //targetStack_1.RemoveTarget(&t1);
    //targetStack_1.ListStack();    

    //targetStack_2.RemoveTarget(&ts_7);
    //targetStack_2.RemoveTarget(&ts_1);
    //targetStack_2.ListStack();

    //TVector3 n = targetStack_1.GetNormalVector_ofStack();
    //printf("n.X() = %f\tn.Y() = %f\tn.Z() = %f\n", n.X(), n.Y(), n.Z());    
    //printf("n.R() = %f\tn.Theta() = %f\tn.Phi() = %f\n", n.Mag(), TMath::RadToDeg()*n.Theta(), TMath::RadToDeg()*n.Phi());    
    targetStack_1.SetStackTiltX(30);
    targetStack_1.SetStackTiltY(45);
    //n = targetStack_1.GetNormalVector_ofStack();
    //printf("\nOrientation\n");
    //printf("n.X() = %f\tn.Y() = %f\tn.Z() = %f\n", n.X(), n.Y(), n.Z());    
    //printf("n.R() = %f\tn.Theta() = %f\tn.Phi() = %f\n", n.Mag(), TMath::RadToDeg()*n.Theta(), TMath::RadToDeg()*n.Phi());    

    //______________________Define beam and particles as UNISIon_______________________________
    //9Li beam, 75 MeV,
    int Zb = 3, Ab = 9;
    double kinE_beam = 75;
    double mass_beam = gNucData->get_mass_Z_A(Zb, Ab);
    double total_energy = kinE_beam + mass_beam;
    double momentum = sqrt(pow(total_energy, 2) - pow(mass_beam, 2));
    double theta = 0;
    double phi = 0;
    TLorentzVector q(momentum*sin(theta)*cos(phi),momentum*sin(theta)*sin(phi),momentum*cos(theta),total_energy);

    UNISIon fBeam;
    fBeam.fMomentum = q;
    fBeam.fZ = Zb;
    fBeam.fA = Ab;
    fBeam.fMass = mass_beam;

    UNISIon part1, part2;
    int Z1 = 3, A1 = 9, Z2=5, A2 = 11;
    double kinE1 = 64.13, kinE2 = 175.87;
    double mass1 = gNucData->get_mass_Z_A(Z1, A1);
    double mass2 = gNucData->get_mass_Z_A(Z2, A2);
    double totE1 = kinE1 + mass1;
    double totE2 = kinE2 + mass2;
    double mom1 = sqrt(pow(totE1, 2) - pow(mass1, 2)), mom2 = sqrt(pow(totE2, 2) - pow(mass2, 2));
    double theta1 = 92.87, theta2 = 130.5;
    double phi1 = 15, phi2 = 195;
    double th1_rad = TMath::DegToRad()*theta1;
    double ph1_rad = TMath::DegToRad()*phi1;
    double th2_rad = TMath::DegToRad()*theta2;
    double ph2_rad = TMath::DegToRad()*phi2;
    TLorentzVector q1(mom1*sin(th1_rad)*cos(ph1_rad),mom1*sin(th1_rad)*sin(ph1_rad),mom1*cos(th1_rad),totE1);
    TLorentzVector q2(mom2*sin(th2_rad)*cos(ph2_rad),mom2*sin(th2_rad)*sin(ph2_rad),mom2*cos(th2_rad),totE2);

    part1.fMomentum = q1;
    part1.fZ = Z1;
    part1.fA = A1;
    part1.fMass = mass1;
    part2.fMomentum = q2;
    part2.fZ = Z2;
    part2.fA = A2;
    part2.fMass = mass2;

//__________________Check propagation of the beam and particles_________________________________________________________


    for(int i = 0; i < 30; ++i)
    {
        std::vector<UNISIon> particles = {part1, part2};
        //Print status before anything entered the target
        printf("\n================================================\n");
        printf("Beam before entering the target\nZ_beam = %d\tA_beam = %d\nkinE_beam = %f\n", Zb, Ab, kinE_beam);
        
        Interaction_Point p1 = targetStack_1.GetInteractionPoint();
        std::setprecision(20);
        printf("p1.target_number = %d\np1.depth = %f\n", p1.target_number, p1.depth);

        Interaction_Point p2 = targetStack_2.GetInteractionPoint();
        printf("\np2.target_number = %d\np2.depth = %f\n", p2.target_number, p2.depth);

        UNISIon fBeam_preInt_1 = targetStack_1.GetBeam_preInteraction(fBeam, p1);

        UNISIon fBeam_preInt_2 = targetStack_2.GetBeam_preInteraction(fBeam, p2);

        UNISIon fBeam_test = testStack.GetBeam_preInteraction(fBeam, p0);
        printf("Target thickness in um = %f\n", testStack.GetTarget_fromStack(p0.target_number).GetTargetThickness_um());
        printf("Beam after degradation in testStack, before interaction:\nZ_beam = %d\tA_beam = %d\nkinE_beam = %f\ntheta_angle = %f\n\n",
                fBeam_test.fZ, fBeam_test.fA, fBeam_test.fMomentum.E() - fBeam_test.fMomentum.M(), TMath::RadToDeg()*fBeam_test.fMomentum.Vect().Theta());

        //Print status after beam degradation, before interaction
        printf("Target thickness in um = %f\n", targetStack_1.GetTarget_fromStack(p1.target_number).GetTargetThickness_um());
        printf("Beam after degradation in target stack 1, before interaction:\nZ_beam = %d\tA_beam = %d\nkinE_beam = %f\ntheta_angle = %f\n\n",
                fBeam_preInt_1.fZ, fBeam_preInt_1.fA, fBeam_preInt_1.fMomentum.E() - fBeam_preInt_1.fMomentum.M(), TMath::RadToDeg()*fBeam_preInt_1.fMomentum.Vect().Theta());

        //Print status after beam degradation, before interaction
        printf("Beam after degradation in target stack 2, before interaction:\nZ_beam = %d\tA_beam = %d\nkinE_beam = %f\n\n",
                fBeam_preInt_2.fZ, fBeam_preInt_2.fA, fBeam_preInt_2.fMomentum.E() - fBeam_preInt_2.fMomentum.M());

        printf("\n================================================\n");
        printf("Particles created in interaction, before target degradation\n");
        printf("Z_part1 = %d\tA_part1 = %d\nkinE_part1 = %f\ntheta1=%f\n", Z1, A1, kinE1, theta1);
        printf("\nZ_part2 = %d\tA_part2 = %d\nkinE_part2 = %f\ntheta2=%f\n\n", Z2, A2, kinE2, theta2);

        //double eff_Thickness_part1 = fabs((targetStack_1.GetTarget_fromStack(p1.target_number).GetTargetThickness_um() - p1.depth)/cos(th1_rad));
        //double eff_Thickness_part2 = fabs((targetStack_1.GetTarget_fromStack(p1.target_number).GetTargetThickness_um() - p1.depth)/cos(th2_rad));
        double eff_Thickness_part1 = fabs((p1.depth)/cos(th1_rad));
        double eff_Thickness_part2 = fabs((p1.depth)/cos(th2_rad));
        printf("Eff_thickness_part1: %f\n", eff_Thickness_part1);
        printf("Eff_thickness_part2: %f\n", eff_Thickness_part2);

        
        std::vector<UNISIon> parts_afterDegrad = targetStack_1.PropagateParticles(particles, p1);
        printf("=============== STACK 1 ==================\n");
        printf("Particles created in interaction, after target degradation\n");
//
        for(auto& part : parts_afterDegrad)
        {
            printf("\nZ = %d\tA = %d\nkinE = %f\n", part.fZ, part.fA, part.fMomentum.E()-part.fMomentum.M());
        }

        parts_afterDegrad = targetStack_2.PropagateParticles(particles, p2);
        printf("=============== STACK 2 ==================\n");
        printf("Particles created in interaction, after target degradation\n");

        for(auto& part : parts_afterDegrad)
        {
            printf("\nZ = %d\tA = %d\nkinE = %f\n", part.fZ, part.fA, part.fMomentum.E()-part.fMomentum.M());
        }

        getchar();
    }

    return 0;
}