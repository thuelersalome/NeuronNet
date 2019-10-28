#include "network.h"
#include "random.h"


void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}
std::pair<size_t, double> Network::degree(const size_t& n) const{
	std::vector<std::pair<size_t, double>> neighbors_neurons(neighbors(n));
	double summ;
	for(auto neighbor : neighbors_neurons){
		summ+=neighbor.second;
	}
	return {neighbors_neurons.size(),summ};
}
std::vector< std::pair< size_t, double > > 	Network::neighbors (const size_t & indice) const
{
	std::vector< std::pair< size_t, double > > connected_neurons;
	for(std::map<std::pair<size_t,size_t>,double>::const_iterator i=links.lower_bound({indice,0}); i!=links.end() and ((i->first).first)==indice; ++i){
		std::pair<size_t, double> insert_neurons ((i->first).second, i->second);
		connected_neurons.push_back(insert_neurons);
		/*
		if(i->first.first==indice)
		{
			connected_neurons.push_back({i->first.second,i->second});
		}
		else if(i->first.second==indice)
		{
			connected_neurons.push_back({i->first.first,i->second});
		}*/
	}
	return connected_neurons;
}

std::set<size_t> Network::step(const std::vector<double>& thalamus_intens)
{
	std::set<size_t> firing_neurons;
	for(size_t i(0);i<neurons.size();++i){ //i est un neuron
		double intensity(0);
		if(neurons[i].is_inhibitory()){
			intensity+=0.4*thalamus_intens[i];
		}
		else{
			intensity+=thalamus_intens[i];
		}
		std::vector< std::pair< size_t, double > > linked_neurons (neighbors(i));
		for(size_t n(0); n<linked_neurons.size();++n){ //n est un neuron connecté (neighbor)
			if (neurons[linked_neurons[n].first].firing()){
				firing_neurons.insert(n);
				intensity+=linked_neurons[n].second;
			}
		}
		neurons[i].input(intensity);
		if(neurons[i].firing()){
			
			neurons[i].reset();
		}
		else{
			
			neurons[i].step();
		}
			
	}
	return firing_neurons;		
}
