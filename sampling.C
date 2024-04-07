#include <iostream>
#include <TVector2.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>

using namespace std;

double diskToSphereRatio = 10.0;

#include <TRandom3.h>

TRandom3 *rnd = new TRandom3(0);

TVector2 PointOnDisk() {
    double r = TMath::Sqrt(rnd->Uniform());
    double theta = rnd->Uniform() * 2 * TMath::Pi();
    return TVector2(r * TMath::Cos(theta), r * TMath::Sin(theta));
}

pair<bool, TVector3>
IntersectionLineSphere(const TVector3 &lineOrigin, const TVector3 &lineDirection) {
    // sphere origin is always (0,0)
    // return the first intersection point
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    const double a = lineDirection.Dot(lineDirection);
    const double b = 2 * lineDirection.Dot(lineOrigin);
    const double c = lineOrigin.Dot(lineOrigin) - 1;

    const double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {false, {0, 0, 0}};
    }
    const double t1 = (-b + sqrt(discriminant)) / (2 * a);
    const double t2 = (-b - sqrt(discriminant)) / (2 * a);
    const double t = min(t1, t2);
    return {true, lineOrigin + t * lineDirection};
}

TF1 zenithAngle("zenithAngle", "TMath::Power(TMath::Cos(x), 2)", 0, TMath::Pi() / 2.0);
TF1 zenithAngleCorrected("zenithAngleCorrected", "zenithAngle(x) / TMath::Cos(x)", 0, TMath::Pi() / 2.0);

void sampling() {
    TH1D *zenithHist = new TH1D("Zenith", "Zenith", 100, 0, TMath::Pi() / 2.0);
    TH1D *zenithCompareHist = new TH1D("Zenith Compare", "Zenith Compare", 100, 0, TMath::Pi() / 2.0);
    TH2D *positionInDisk = new TH2D("Position in Disk", "Position in Disk", 1000, -10, 10, 1000, -10, 10);
    TH1D *distanceToCenter = new TH1D("distanceToCenter", "distanceToCenter", 1000, 0, 5);
    while (zenithHist->GetEntries() < 10000000) {
        break;
        const TVector2 pointOnDisk = PointOnDisk() * diskToSphereRatio;
        const TVector3 point(pointOnDisk.X(), 1, pointOnDisk.Y());

        const double phi = rnd->Uniform() * TMath::TwoPi();
        const double zenith = zenithAngle.GetRandom();

        const TVector3 direction = {TMath::Sin(zenith) * TMath::Cos(phi), -1.0 * TMath::Cos(zenith),
                                    TMath::Sin(zenith) * TMath::Sin(phi)};


        const auto [intersectionFlag, intersection] = IntersectionLineSphere(point, direction);
        if (!intersectionFlag) {
            continue;
        }
        positionInDisk->Fill(pointOnDisk.X(), pointOnDisk.Y());
        distanceToCenter->Fill(pointOnDisk.Mod());
        zenithHist->Fill(zenith);
        zenithCompareHist->Fill(zenithAngleCorrected.GetRandom());
    }

    while (zenithHist->GetEntries() < 10000000) {
        // break;
        const auto zenith = zenithAngleCorrected.GetRandom();
        // ellipse center is (tan(zenith), 0) and major axis is 1/cos(zenith), minor axis is 1
        // random point in ellipse, uniformly distributed

        // TODO: improve with smarter sampling
        double y = 1;
        const double xOffset = TMath::Tan(zenith);
        double x = 1 / TMath::Cos(zenith) + xOffset;
        while (TMath::Power(x - TMath::Sin(zenith), 2) * TMath::Cos(zenith) + TMath::Power(y, 2) > 1) {
            y = rnd->Uniform(-1, 1);
            x = rnd->Uniform(-1, 1) * 1 / TMath::Cos(zenith) + xOffset;
        }
        const double phi = rnd->Uniform() * TMath::TwoPi();
        TVector2 positionOnDisk(x, y);
        // rotate by phi
        x = positionOnDisk.X() * TMath::Cos(phi) - positionOnDisk.Y() * TMath::Sin(phi);
        y = positionOnDisk.X() * TMath::Sin(phi) + positionOnDisk.Y() * TMath::Cos(phi);
        positionOnDisk = TVector2(x, y);
        positionInDisk->Fill(positionOnDisk.X(), positionOnDisk.Y());
        distanceToCenter->Fill(positionOnDisk.Mod());
        zenithHist->Fill(zenith);
    }

    // normalize histogram
    // zenithHist->Scale(1.0 / zenithHist->Integral());

    // canvas
    TCanvas *c = new TCanvas("c", "c", 800, 600);

    distanceToCenter->Draw();
    // positionInDisk->Draw("colz");

    // zenithHist->SetLineColor(kRed);
    // zenithHist->Draw();


    // draw zenithAngle and zenithAngleCorrected on top
    // zenithAngle.SetLineColor(kGreen);
    // zenithCompareHist->Draw("same");

}
