

void plot()
{
	ifstream in("data.dat");
	int id;
	double x, y;

	TH1D *h = new TH1D("h", "h", 100, 0.0, 50.0);
	while(in >> id >> x >> y) {
		double vel = sqrt(x*x + y*y);
		h->Fill(vel);
	}
	in.close();
	h->Draw();
	}

