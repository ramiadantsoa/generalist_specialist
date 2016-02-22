#ifndef __PATCH_H_
#define __PATCH_H_
class Patch {
    public:
        Patch(double _x, double _y, double _q, int _M) {
            this->x = _x;
            this->y = _y;
            this->quality = _q;
            this->M = _M; // The number of species, i.e., length of occupancy matrix
            this->occupancy = new int[_M];
            this->colonization = new double[_M];
			this->fitness = new double[_M];
            for(int i = 0; i<_M; i++) {
                this->occupancy[i] = 0;
                this->colonization[i] = 0.0;
				this->fitness[i]= i ==0 ? 1 : exp((_M-1)*(_M-1)/2.0*cos(2.0*PI*(_q - (double)i /(double)(_M-1)))- 6.058104255427814);
            }
        }

        ~Patch() {
            delete[] this->occupancy;
            delete[] this->colonization;
			delete[] this->fitness;
        }

		double get_fitness(int sp_id){
			double fitn = sp_id==0 ? 1 : this->fitness[sp_id];
			return fitn;
		}

        double get_x() { 
            return this->x; 
        }

        double get_y() { 
            return this->y; 
        }

        double get_quality() {
            return this->quality;
        }

        int count_present_species() {
            int count = 0;
            for(int i = 0; i < this->M; i++) {
                if (this->occupancy[i] == 1)
                    count += 1;
            }
            return count;
        }

		double get_colonization(int m){
			return this->colonization[m];
		}

		void modify_col_patch(int m, double value){
			this->colonization[m]+=value;
			if( this->colonization[m]<0.0 || this->colonization[m]<chop_0){
				this->colonization[m]= 0.0;
			}
		}

		void set_col_to_0(int m){
			this->colonization[m]=0.0;
		}

        int *occupancy;
    private:
        int M;
        double x, y; 
        double quality; 
		double *fitness;
		double *colonization;
};

void output_occupancy(Patch *p, int M){
	for( int m = 0; m<M ; m++){
		cout<<p->occupancy[m]<<" \t ";
	}
	cout<<endl;
}

#endif // __PATCH_H_