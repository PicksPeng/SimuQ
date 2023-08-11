from simuq.solver import generate_as

def to_bin(res, n) :
    b = bin(res)[2:]
    return '0' * (n - len(b)) + b

class BaseProvider() :
    def __init__(self) :
        self.prog = None
        self.task = None
        self.qs_names = None

    def print_sites(self) :
        if self.prog == None :
            raise Exception("No compiled job in record.")
        print("Order of sites:", self.qs_names)

class BraketProvider(BaseProvider) :
    def __init__(self) :

        try :
            from braket.aws import AwsDevice, AwsQuantumTask
        except :
            raise Exception("Please install amazon braket python sdk: 'pip install amazon-braket-sdk' \n and configure your region and credentials according to instructions in: https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html#configuration")

        self.backend_aais = dict()
        self.backend_aais[('ionq', 'harmony')] = ['heis_aais']
        self.backend_aais[('ionq', 'aria-1')] = ['heis_aais']
        self.backend_aais[('quera', 'Aquila')] = ['rydberg1d_global', 'rydberg2d_global']

        super().__init__()

    def supported_backends(self) :
        for (comp, dev) in self.backend_aais.keys() :
            print(f'Hardware provider: {comp}  -  Device name: {dev}  -  AAIS supports: {self.backend_aais_supports[(comp, dev)]}')

    def compile(self, qs, provider, device, aais, tol = 0.01, trotter_num = 6, verbose = 0) :
        if (provider, device) not in self.backend_aais.keys() :
            raise Exception("Not supported hardware provider or device.")
        if aais not in self.backend_aais[(provider, device)] :
            raise Exception("Not supported AAIS on this device.")

        self.qs_names = qs.print_sites()

        self.provider = provider
        self.device = device

        if self.provider == 'quera' :
            
            nsite = qs.num_sites
            
            if aais == 'rydberg1d_global' :
                from simuq.aais.rydberg1d_global import Rydberg1DGlobalQMachineFactory
                from simuq.backends.braket_rydberg1d_global import transpile
                mach = Rydberg1DGlobalQMachineFactory.generate_qmachine(nsite)
                comp = transpile
            elif aais == 'rydberg2d_global' :
                from simuq.aais.rydberg2d_global import Rydberg2DGlobalQMachineFactory
                from simuq.backends.braket_rydberg2d_global import transpile
                mach = Rydberg2DGlobalQMachineFactory.generate_qmachine(nsite)
                comp = transpile

            layout, sol_gvars, boxes, edges = generate_as(qs, mach, trotter_num = 1,
                                                          solver = 'least_squares',
                                                          solver_args = {'tol' : tol, 'with_time_penalty' : True},
                                                          override_layout = [i for i in range(nsite)],
                                                          verbose = verbose)
            
            self.ahs_prog = comp(sol_gvars, boxes, edges, verbose = verbose)
            self.prog = self.ahs_prog
            self.layout = layout

    def visualize_quera(self) :
        from braket.ahs.atom_arrangement import AtomArrangement, SiteType
        import matplotlib.pyplot as plt
        import numpy as np

        braket_prog = self.ahs_prog

        def show_register(register):
            filled_sites = [site.coordinate for site in register._sites if site.site_type == SiteType.FILLED]
            empty_sites = [site.coordinate for site in register._sites if site.site_type == SiteType.VACANT]
            fig = plt.figure(figsize=(5, 5))
            if len(filled_sites) > 0:
                plt.plot(np.array(filled_sites)[:, 0], np.array(filled_sites)[:, 1], 'r.', ms=15, label='filled')
            if len(empty_sites) > 0:
                plt.plot(np.array(empty_sites)[:, 0], np.array(empty_sites)[:, 1], '.', color='k', ms=2, label='vacant')
            plt.legend(bbox_to_anchor=(1.1, 1.05))

        show_register(braket_prog.register)

        def show_global_drive(drive):
            data = {
                'amplitude [rad/s]': drive.amplitude.time_series,
                'detuning [rad/s]': drive.detuning.time_series,
                'phase [rad]': drive.phase.time_series,
            }
    
            fig, axes = plt.subplots(3, 1, figsize=(5, 4), sharex=True)
            for ax, data_name in zip(axes, data.keys()):
                if data_name == 'phase [rad]':
                    ax.step(data[data_name].times(), data[data_name].values(), '.-', where='post')
                else:
                    ax.plot(data[data_name].times(), data[data_name].values(), '.-')
                ax.set_ylabel(data_name)
                ax.grid(ls=':')
            axes[-1].set_xlabel('time [s]')
            plt.tight_layout()

        show_global_drive(braket_prog.hamiltonian)
        plt.show()

    def visualize(self) :
        if self.prog == None :
            raise Exception("No compiled job in record.")
        if self.provider == 'quera' :
            self.visualize_quera()

    def run(self, shots = 1000, on_simulator = False) :
        if self.prog == None :
            raise Exception("No compiled job in record.")
        
        if self.provider == 'quera' :
            from braket.aws import AwsDevice, AwsQuantumTask
            from braket.devices import LocalSimulator
            if on_simulator :
                simulator = LocalSimulator("braket_ahs")
                self.task = simulator.run(self.ahs_prog, shots = shots)
                meta = self.task.metadata()
                print("Submitted.")
            else :
                aquila_qpu = AwsDevice("arn:aws:braket:us-east-1::device/qpu/quera/" + device)
                self.task = aquila_qpu.run(self.ahs_prog, shots = shots)
                print("Submitted.")
                meta = self.task.metadata()
                print("Task arn: ", meta['quantumTaskArn'])
                print("Task status: ", meta['status'])

    def results(self, task_arn = None) :
        if task_arn != None :
            from braket.aws import AwsQuantumTask
            task = AwsQuantumTask(arn=task_arn)
        else :
            if self.task == None :
                raise Exception("No submitted job in record.")
            task = self.task

        state = task.state()
        if state != 'COMPLETED' :
            print("Job is not completed")
            print(task.metadata())
            return

        result = task.result()

        N = len(self.layout)
        counts = dict([])
        for i in range(result.task_metadata.shots) :
            if sum(result.measurements[i].pre_sequence) == N :
                s = ''
            for j in range(N) :
                s += 'g' if result.measurements[i].post_sequence[j] == 1 else 'r'
            if s not in counts.keys() :
                counts[s] = 1
            else :
                counts[s] +=1

        def int_to_gr(N, s) :
            ret = ''
            for i in  range(N) :
                ret += 'g' if (s>>i)&1 == 0 else 'r'
            return ret

        def extract_prob(counts) :
            tot = sum(counts.values())
            freq = dict([])
            for k in counts.keys() :
                freq[k] = counts[k] / tot
            exp_prob = []
            for i in range(1<<N) :
                key = int_to_gr(N, i)
                if key not in freq.keys() :
                    exp_prob.append(0)
                else :
                    exp_prob.append(freq[key])
            return exp_prob

        prob = extract_prob(counts)
        ret = dict()
        for i in range(1<<N) :
            if abs(prob[i]) > 1e-6 :
                ret[to_bin(i, N)] = prob[i]
        return ret





class IonQProvider(BaseProvider) : 
    def __init__(self, API_key = None, from_file = None) :
        if API_key == None :
            if from_file == None :
                raise Exception("No API_key provided.")
            else :
                with open(from_file, 'r') as f :
                    API_key = f.readline()
        self.API_key = API_key
        self.all_backends = ['harmony', 'aria-1', 'aria-2', 'forte']
        
        super().__init__()

    def supported_backends(self) :
        print(self.all_backends)

    def compile(self, qs, backend = 'aria-1', aais = 'heis_aais', tol = 0.01, trotter_num = 6, verbose = 0) :
        
        if backend == 'harmony' :
            nsite = 11
        elif backend == 'aria-1' :
            nsite = 23
        elif backend == 'aria-2' :
            nsite = 23
        elif backend == 'forte' :
            nsite = 32
        elif isinstance(backend, int) :
            if backend > 0 :
                nsite = backend
        else :
            raise Exception("Backend is not supported.")

        if qs.num_sites > nsite :
            raise Exception("Device has less sites than the target quantum system.")
        
        if aais == 'heis_aais' :

            from simuq.aais.heisenberg import HeisenbergQMachineFactory
            from simuq.backends.ionq import transpile

            mach = HeisenbergQMachineFactory.generate_qmachine(qs.num_sites, E = None)
            comp = transpile

        layout, sol_gvars, boxes, edges = generate_as(qs, mach, trotter_num,
                                                      solver = 'least_squares', solver_args = {'tol':tol},
                                                      override_layout = [i for i in range(qs.num_sites)],
                                                      verbose = verbose)
        self.prog = comp(qs.num_sites, sol_gvars, boxes, edges, backend = 'qpu.' + backend, noise_model = backend)
        self.layout = layout
        self.qs_names = qs.print_sites()


    def run(self, shots = 4096, on_simulator = False) :

        import requests, json

        headers = {
            'Authorization': 'apiKey ' + self.API_key,
            'Content-Type': 'application/json',
        }

        data = self.prog.copy()

        data['shots'] = shots
        if on_simulator :
            data['target'] = 'simulator'
        
        #print(data)

        response = requests.post('https://api.ionq.co/v0.3/jobs', headers=headers, data=json.dumps(data))
        self.task = response.json()

        print(self.task)

    def results(self, job_id = None) :
        if job_id == None :
            if self.task != None :
                job_id = self.task['id']
            else :
                raise Exception("No submitted job in record.")

        import requests, json

        headers = {
            'Authorization': 'apiKey ' + self.API_key,
        }

        response = requests.get('https://api.ionq.co/v0.2/jobs/' + job_id, headers=headers)
        res = response.json()

        if res['status'] != 'completed' :
            print("Job is not completed")
            print(res)
            return

        def layout_rev(res) :
            n = len(self.layout)
            b = to_bin(res, n)
            ret = ''
            for i in range(n) :
                ret += b[self.layout[i]]
            return ret

        def results_from_data(data) :
            ret = dict()
            for key in data.keys() :
                ret[layout_rev(int(key))] = data[key]
            return ret

        return results_from_data(res['data']['histogram'])
        


class QuTiPProvider(BaseProvider) :
    def __init__(self) :
        try :
            import qutip as qp
        except :
            raise Exception("No QuTiP package detected. Please try command 'pip install simuq[qutip]'. ")
        self.prog = None
        self.fin = None

    def compile(self, qs) :
        import qutip as qp
        self.n = qs.num_sites
        self.init = qp.basis(1<<self.n)
        self.prog = (qs.to_qutip(), qs.total_time())
        self.qs_names = qs.print_sites()
        print("Compiled.")
        #return self.prog
        
    def run(self, shots = None, on_simulator = None) :
        import qutip as qp
        if self.prog == None :
            raise Exception("No compiled job in record.")
        self.fin = qp.sesolve(self.prog[0], self.init, [0, self.prog[1]])
        print("Solved.")
        #return self.fin

    def results(self) :
        import numpy as np
        if self.fin == None :
            raise Exception("No submitted job in record.")
        self.res = dict()
        for i in range(1<<self.n) :
            self.res[to_bin(i, self.n)] = np.abs(self.fin.states[1][i][0][0]) ** 2
        return self.res
