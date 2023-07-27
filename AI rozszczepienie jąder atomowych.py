import numpy as np
import matplotlib.pyplot as plt
import random
import unittest
from weapon import Weapon
from physics import Power, Momentum, AngularMomentum, Torque, MomentOfInertia


class FissionAnalysis:
    def __init__(self):
        self.neutron_mass = 1.008664916  # masa atomowa neutronu w u
        self.daughter_mass = 0.92828918  # masa atomowa jądra po rozszczepieniu w u

    def fission(self, mass):
        """
        Funkcja rozszczepienia jądra atomowego.

        Argumenty:
            mass (float): Masa jądra atomowego w u.

        Zwraca:
            float: Energia rozszczepienia jądra atomowego w MeV.
        """
        energy = (mass - self.neutron_mass - self.daughter_mass) * 931.5  # energia rozszczepienia w MeV

        # Oblicz liczbę neutronów emitowanych podczas rozszczepienia.
        neutron_count = (mass - self.neutron_mass - self.daughter_mass) / self.neutron_mass

        # Oblicz liczbę produktów rozszczepienia.
        product_count = 2 * neutron_count + 1

        # Wypisz informacje o rozszczepieniu.
        print(f"Masa jądra: {mass} u")
        print(f"Energia rozszczepienia: {energy} MeV")
        print(f"Liczba neutronów emitowanych: {neutron_count}")
        print(f"Liczba produktów rozszczepienia: {product_count}")

        # Sprawdź, czy rozszczepienie jest możliwe.
        if mass < self.critical_mass:
            print("Rozszczepienie nie jest możliwe.")
            return None

        # Wygeneruj produkty rozszczepienia.
        products = []
        for i in range(product_count):
            product = random.choice(self.products)
            products.append(product)

        # Wypisz produkty rozszczepienia.
        print("Produkty rozszczepienia:")
        for product in products:
            print(product)

        return energy

    def collect_data(self, start_mass, end_mass, step=1):
        """
        Funkcja zbiera dane z rozszczepienia jądra atomowego.

        Argumenty:
            start_mass (float): Masa początkowa jądra atomowego w u.
            end_mass (float): Masa końcowa jądra atomowego w u.
            step (float): Krok między masami w u.

        Zwraca:
            (list, list): Lista mas i lista energii.
        """
        masses = np.arange(start_mass, end_mass + 1, step)
        energies = []
        for mass in masses:
            energy = self.fission(mass)
            if energy is not None:
                energies.append(energy)
        return masses, energies

        # Zapisz dane do pliku.
    with open("data.csv", "w") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(["Mass", "Energy"])
        for mass, energy in zip(masses, energies):
            writer.writerow([mass, energy])

        # Wykreśl dane.
        plt.plot(masses, energies)
        plt.xlabel("Masa w u")
        plt.ylabel("Energia w MeV")
        plt.show()


    def plot_energy_vs_mass(self, masses, energies):
        """
        Funkcja rysuje wykres zależności energii rozszczepienia od masy jądra atomowego.

        Argumenty:
            masses (list): Lista mas jąder atomowych w u.
            energies (list): Lista energii rozszczepienia w MeV.

        Zwraca:
            matplotlib.figure.Figure: Obiekt wykresu.
        """
        fig, ax = plt.subplots()
        ax.plot(masses, energies, linestyle='dotted', marker='o', color='b')
        ax.set_xlabel("Masa jądra atomowego (u)")
        ax.set_ylabel("Energia rozszczepienia (MeV)")
        ax.set_title("Zależność energii rozszczepienia od masy jądra atomowego")
        ax.grid(True)

        # Dodaj opis do wykresu.
        ax.annotate("Energia rozszczepienia rośnie wraz ze wzrostem masy jądra atomowego.",
                    xy=(0.5, 0.5),
                    xytext=(-10, 10),
                    textcoords="offset points",
                    size="10",
                    ha="center",
                    va="center")

        # Zapisz wykres do pliku.
        plt.savefig("energy_vs_mass.png")

        return fig


    def get_all_isotopes(element):
        """
        Zwraca listę wszystkich izotopów podanego pierwiastka.

        Argumenty:
            element (str): Pierwiastek.

        Zwraca:
            list: Lista izotopów.
        """
        # Pobierz listę wszystkich izotopów z bazy danych.
        isotopes = get_isotopes_from_database(element)

        # Posortuj listę izotopów według masy atomowej.
        isotopes.sort(key=lambda x: x.atomic_mass)

        # Dodaj do listy izotopów informacje o trwałości.
        for isotope in isotopes:
            isotope.half_life = get_half_life_from_database(isotope.id)

        # Filtruj listę izotopów, pozostawiając tylko izotopy stabilne.
        stable_isotopes = []
        for isotope in isotopes:
            if isotope.half_life is not None:
                stable_isotopes.append(isotope)

        return stable_isotopes



    def get_isotope_by_mass(element, mass):
        """
        Zwraca izotop o podanej masie atomowej dla podanego pierwiastka.

        Argumenty:
            element (str): Pierwiastek.
            mass (float): Masa atomowa.

        Zwraca:
            Isotope: Izotop.
        """
        # Pobierz listę wszystkich izotopów dla podanego pierwiastka.
        isotopes = get_all_isotopes(element)

        # Znajdź izotop o podanej masie atomowej.
        isotope = next(filter(lambda x: x.atomic_mass == mass, isotopes), None)

        # Jeśli izotop nie został znaleziony, zgłoś błąd.
        if isotope is None:
            raise ValueError(f"No isotope with mass {mass} was found for element {element}.")

        # Pobierz inne informacje o izotopie, takie jak okres połowicznego rozpadu, energia wiązania jądrowego i ilość neutronów.
        isotope.half_life = get_half_life_from_database(isotope.id)
        isotope.binding_energy = get_binding_energy_from_database(isotope.id)
        isotope.neutron_count = get_neutron_count_from_database(isotope.id)

        return isotope



    def get_fission_energy(isotope):
        """
        Zwraca energię rozszczepienia podanego izotopu.

        Argumenty:
            isotope (Isotope): Izotop.

        Zwraca:
            float: Energia rozszczepienia.
        """
        # Oblicz energię rozszczepienia.
        energy = (
        isotope.atomic_mass
        - isotope.neutron_mass
        - isotope.daughter_mass
    ) * 931.5 

        # Dodaj do energii rozszczepienia energię kinetyczną produktów rozszczepienia.
        for product in isotope.products:
            energy += product.kinetic_energy

        # Oblicz ilość neutronów emitowanych podczas rozszczepienia.
        neutron_count = (
        isotope.atomic_mass
        - isotope.daughter_mass
        - isotope.neutron_mass
    ) / neutron_mass

        # Wypisz informacje o rozszczepieniu.
        print(f"Energia rozszczepienia: {energy} MeV")
        print(f"Liczba neutronów emitowanych podczas rozszczepienia: {neutron_count}")

        return energy


       
    def plot_fission_energy_vs_mass(elements):
        """
        Rysuje wykres zależności energii rozszczepienia od masy atomowej dla podanych pierwiastków.

        Argumenty:
            elements (list): Lista pierwiastków.
        """
        # Pobierz listę wszystkich izotopów dla podanych pierwiastków.
        isotopes = [get_all_isotopes(element) for element in elements]

        # Połącz listy izotopów w jedną listę.
        isotopes = sum(isotopes, [])

        # Posortuj listę izotopów według masy atomowej.
        isotopes.sort(key=lambda x: x.atomic_mass)

        # Utwórz wykres.
        plt.figure()

        # Dodaj do wykresu dane.
        for isotope in isotopes:
            plt.plot(isotope.atomic_mass, isotope.fission_energy, marker='o', color='b')

        # Dodaj etykiety do osi wykresu.
        plt.xlabel("Masa atomowa (u)")
        plt.ylabel("Energia rozszczepienia (MeV)")

        # Dodaj tytuł wykresu.
        plt.title("Zależność energii rozszczepienia od masy atomowej")

        # Pokaż wykres.
        plt.show()

        # Dodaj opis do wykresu.
        plt.annotate("Energia rozszczepienia rośnie wraz ze wzrostem masy atomowej.",
                    xy=(0.5, 0.5),
                    xytext=(-10, 10),
                    textcoords="offset points",
                    size="10",
                    ha="center",
                    va="center")

        # Zapisz wykres do pliku.
        plt.savefig("energy_vs_mass.png")

        return fig

        # Dodaj do wykresu linie trendu.
    plt.plot(isotopes.atomic_mass, isotopes.fission_energy, '-r')

        # Dodaj etykiety do linii trendu.
    plt.annotate("Linia trendu",
                xy=(0.5, 0.5),
                xytext=(-10, 10),
                textcoords="offset points",
                size="10",
                ha="center",
                va="center")

        # Pokaż wykres.
    plt.show()

        # Zapisz wykres do pliku.
    plt.savefig("energy_vs_mass_with_trend_line.png")

class Neuron:
    def __init__(self, weights, biases, activation_function):
        self.weights = weights
        self.biases = biases
        self.activation_function = activation_function

        # Dodaj domyślną funkcję aktywacji.
        if activation_function is None:
            self.activation_function = lambda x: x

    def forward_propagate(self, inputs):
        outputs = np.dot(inputs, self.weights) + self.biases
        outputs = self.activation_function(outputs)

        # Dodaj możliwość zwracania wag i biasów.
        return outputs, self.weights, self.biases

    def backward_propagate(self, inputs, outputs, delta):
        # Oblicz gradient funkcji aktywacji.
        grad_activation = self.activation_function(outputs).derivative()

        # Oblicz gradienty wag i biasów.
        grad_weights = np.dot(inputs.T, grad_activation * delta)
        grad_biases = np.sum(grad_activation * delta, axis=0)

        return grad_weights, grad_biases

    def update_parameters(self, learning_rate, grad_weights, grad_biases):
        # Zaktualizuj wagi i biasy.
        self.weights -= learning_rate * grad_weights
        self.biases -= learning_rate * grad_biases

    def train(self, inputs, outputs, epochs, learning_rate):
        for epoch in range(epochs):
            # Przeprowadź propagację wsteczną dla wszystkich danych treningowych.
            for input, output in zip(inputs, outputs):
                grad_weights, grad_biases = self.backward_propagate(input, output, 1)
                self.update_parameters(learning_rate, grad_weights, grad_biases)

            # Wyznacz dokładność sieci na danych treningowych.
            accuracy = self.evaluate(inputs, outputs)

            # Wypisz dokładność sieci.
            print(f"Epoch {epoch + 1}: {accuracy * 100:.2f}%")

    def evaluate(self, inputs, outputs):
        # Przeprowadź propagację w przód dla wszystkich danych testowych.
        outputs_pred = [self.forward_propagate(input)[0] for input in inputs]

        # Oblicz dokładność sieci.
        accuracy = np.mean(np.argmax(outputs, axis=1) == np.argmax(outputs_pred, axis=1))

        return accuracy



class NeuralNetwork:
    def __init__(self, layers, activation_functions):
        self.layers = layers
        self.activation_functions = activation_functions

    def forward_propagate(self, inputs):
        outputs = inputs
        for layer in self.layers:
            outputs = layer.forward_propagate(outputs)
        return outputs

    def backpropagate(self, inputs, outputs, labels):
        # Calculate the loss
        loss = self.loss_function(outputs, labels)

        # Backpropagate the loss
        gradients = []
        for layer in reversed(self.layers):
            gradients.append(layer.backpropagate(inputs, outputs, labels))

        # Update the weights and biases
        for layer, gradient in zip(self.layers, gradients):
            layer.update_weights(gradient)

    def train(self, inputs, labels, epochs, learning_rate):
        for epoch in range(epochs):
            loss = 0
            for i in range(len(inputs)):
                outputs = self.forward_propagate(inputs[i])
                loss += self.loss_function(outputs, labels[i])

            self.backpropagate(inputs, outputs, labels)

        return loss

    def predict(self, inputs):
        outputs = self.forward_propagate(inputs)
        return outputs
   
class FissionAnalysis:
    def __init__(self, data):
        self.data = data

        # Dodaj domyślne parametry.
        self.dropna_thresh = 0.1
        self.fillna_val = 0.0
        self.normalization_norm = 'l2'

        # Dodaj możliwość obsługi różnych typów danych.
        if not isinstance(data, (pd.DataFrame, np.ndarray)):
            raise ValueError('data musi być typu pandas.DataFrame lub numpy.ndarray')

    def process_data(self):
        """
        Przetwórz dane, aby przygotować je do analizy przez algorytmy uczenia maszynowego.

        Argumenty:
            self (FissionAnalysis): Obiekt klasy FissionAnalysis.

        Zwraca:
            np.ndarray: Przetworzone dane.
        """

        # Usuń puste wiersze.
        self.data = self.data.dropna(thresh=self.dropna_thresh)

        # Zastąp wartości NaN średnią wartościami w kolumnie.
        self.data = self.data.fillna(self.fillna_val)

        # Znormalizuj dane.
        self.data = self.data.astype(float)
        self.data = preprocessing.normalize(self.data, norm=self.normalization_norm)

        return self.data

    
    
    def predict_fission(self, input_data):
        """
        Wykorzystaj model do przewidywania energii rozszczepienia jądra atomowego.

        Argumenty:
            self (FissionAnalysis): Obiekt klasy FissionAnalysis.
            input_data (np.ndarray): Dane wejściowe.

        Zwraca:
            float: Przewidywana energia rozszczepienia jądra atomowego.
        """

        # Przekonwertuj dane wejściowe na wektor.
        input_data = np.array(input_data)

        # Znormalizuj dane wejściowe.
        input_data = preprocessing.normalize(input_data, norm='l2')

        # Przewidz energię rozszczepienia jądra atomowego.
        predicted_fission_energy = model.predict(input_data)

        # Wyznacz dokładność modelu.
        accuracy = np.mean(np.argmax(input_data, axis=1) == np.argmax(predicted_fission_energy, axis=1))

        # Wypisz dokładność modelu.
        print(f"Dokładność: {accuracy * 100:.2f}%")

        # Wyznacz macierz konfuzji.
        confusion_matrix = np.zeros((input_data.shape[1], input_data.shape[1]))
        for i in range(input_data.shape[0]):
            confusion_matrix[np.argmax(input_data[i]), np.argmax(predicted_fission_energy[i])] += 1

        # Wypisz macierz konfuzji.
        print(confusion_matrix)

        return predicted_fission_energy




class FissionModel:
    def __init__(self, nucleus_structure, neutron_energy, number_of_neutrons, temperature, pressure):
        self.nucleus_structure = nucleus_structure
        self.neutron_energy = neutron_energy
        self.number_of_neutrons = number_of_neutrons
        self.temperature = temperature
        self.pressure = pressure

    def calculate_probability_of_fission(self):
        """
        Oblicza prawdopodobieństwo rozszczepienia jądra atomowego.

        Zakładamy, że prawdopodobieństwo rozszczepienia zależy od liczby neutronów,
        temperatury i ciśnienia.

        Zwraca:
            float: Prawdopodobieństwo rozszczepienia jądra atomowego (od 0.0 do 1.0).
        """
        # Zakładamy, że prawdopodobieństwo rozszczepienia zależy od liczby neutronów,
        # temperatury i ciśnienia. Możemy tu zastosować bardziej złożony model,
        # ale w tym przykładzie użyjemy jedynie prostego wzoru.
        probability = (self.number_of_neutrons / 1000 + self.temperature / 500 + self.pressure / 200) / 3

        # Ograniczamy prawdopodobieństwo do przedziału [0.0, 1.0].
        probability = max(0.0, min(1.0, probability))
        return probability

    def calculate_number_of_fission_neutrons(self):
        """
        Oblicza liczbę neutronów powstających w wyniku rozszczepienia jądra atomowego.

        Zakładamy, że liczba neutronów zależy od energii neutronów, liczby neutronów
        w jądrze oraz temperatury.

        Zwraca:
            int: Liczba neutronów powstających w wyniku rozszczepienia.
        """
        # Zakładamy, że liczba neutronów jest proporcjonalna do energii neutronów,
        # liczby neutronów w jądrze i maleje wraz ze wzrostem temperatury.
        number_of_neutrons = int(self.neutron_energy / 10 + self.number_of_neutrons / 50 - self.temperature / 200)

        # Liczba neutronów nie może być ujemna.
        number_of_neutrons = max(0, number_of_neutrons)
        return number_of_neutrons

    def calculate_fission_energy(self, probability_of_fission, number_of_fission_neutrons):
        """
        Oblicza energię rozszczepienia jądra atomowego.

        Argumenty:
            probability_of_fission (float): Prawdopodobieństwo rozszczepienia jądra atomowego (od 0.0 do 1.0).
            number_of_fission_neutrons (int): Liczba neutronów powstałych w wyniku rozszczepienia.

        Zwraca:
            float: Energia rozszczepienia jądra atomowego.
        """
        if not (0.0 <= probability_of_fission <= 1.0):
            raise ValueError("Prawdopodobieństwo rozszczepienia powinno być wartością z zakresu od 0.0 do 1.0.")
        
        if number_of_fission_neutrons <= 0:
            raise ValueError("Liczba neutronów w wyniku rozszczepienia powinna być dodatnia.")

        # Obliczenie masy rozszczepionego jądra atomowego w u
        mass_of_fissioned_nucleus = self.nucleus_structure["mass"] + self.number_of_neutrons * self.nucleus_structure["neutron_mass"]

        # Obliczenie energii rozszczepienia w MeV
        energy = mass_of_fissioned_nucleus * 931.5 * probability_of_fission

        return energy
class FissionControlModel:

    def __init__(self, nucleus_structure, neutron_energy, number_of_neutrons, temperature, pressure):
        self.nucleus_structure = nucleus_structure
        self.neutron_energy = neutron_energy
        self.number_of_neutrons = number_of_neutrons
        self.temperature = temperature
        self.pressure = pressure
        self.neutron_mass = 1.008664916  # masa atomowa neutronu w u
        self.daughter_mass = 0.92828918  # masa atomowa jądra po rozszczepieniu w u


    def get_fission_energy(self):
        """
        Oblicza energię rozszczepienia jądra atomowego.

        Argumenty:
            self (FissionControlModel): Obiekt klasy FissionControlModel.

        Zwraca:
            float: Energia rozszczepienia jądra atomowego.
        """

        # Oblicz prawdopodobieństwo rozszczepienia jądra atomowego.
        probability_of_fission = self.calculate_probability_of_fission()

        # Oblicz liczbę neutronów powstających w wyniku rozszczepienia.
        number_of_fission_neutrons = self.calculate_number_of_fission_neutrons()

        # Oblicz energię rozszczepienia jądra atomowego.
        fission_energy = self.calculate_fission_energy(probability_of_fission, number_of_fission_neutrons)

        return fission_energy

    def calculate_probability_of_fission(self):
        """
        Oblicza prawdopodobieństwo rozszczepienia jądra atomowego.

        Argumenty:
            self (FissionControlModel): Obiekt klasy FissionControlModel.

        Zwraca:
            float: Prawdopodobieństwo rozszczepienia jądra atomowego.
        """

        # Oblicz wartość funkcji penetracji.
        penetration_function_value = self.calculate_penetration_function_value()

        # Oblicz wartość funkcji rozszczepienia.
        fission_function_value = self.calculate_fission_function_value()

        # Oblicz prawdopodobieństwo rozszczepienia jądra atomowego.
        probability_of_fission = penetration_function_value * fission_function_value

        return probability_of_fission

    def calculate_number_of_fission_neutrons(self):
        """
        Oblicza liczbę neutronów powstających w wyniku rozszczepienia.

        Argumenty:
            self (FissionControlModel): Obiekt klasy FissionControlModel.

        Zwraca:
            int: Liczba neutronów powstających w wyniku rozszczepienia.
        """

        # Oblicz wartość funkcji rozszczepienia.
        fission_function_value = self.calculate_fission_function_value()

        # Oblicz liczbę neutronów powstających w wyniku rozszczepienia.
        number_of_fission_neutrons = fission_function_value

        return number_of_fission_neutrons

    def calculate_fission_energy(self, probability_of_fission, number_of_fission_neutrons):
        """
        Oblicza energię rozszczepienia jądra atomowego.

        Argumenty:
            probability_of_fission (float): Prawdopodobieństwo rozszczepienia jądra atomowego (od 0.0 do 1.0).
            number_of_fission_neutrons (int): Liczba neutronów powstałych w wyniku rozszczepienia.

        Zwraca:
            float: Energia rozszczepienia jądra atomowego w MeV.
        """
        if not (0.0 <= probability_of_fission <= 1.0):
            raise ValueError("Prawdopodobieństwo rozszczepienia powinno być wartością z zakresu od 0.0 do 1.0.")
        
        if number_of_fission_neutrons <= 0:
            raise ValueError("Liczba neutronów w wyniku rozszczepienia powinna być dodatnia.")

        # Obliczenie masy rozszczepionego jądra atomowego w u
        mass_of_fissioned_nucleus = self.daughter_mass + self.neutron_mass * number_of_fission_neutrons

        # Obliczenie energii rozszczepienia w MeV
        energy = mass_of_fissioned_nucleus * 931.5 * probability_of_fission

        return energy
class NeutronTracker:

    def __init__(self, reactor):
        self.reactor = reactor

    def track_neutron(self, neutron):
        """
        Śledzi neutron w reaktorze.

        Argumenty:
            neutron (Neutron): Neutron do śledzenia.

        Zwraca:
            float: Energię neutronu po śladzie.
        """

        # Oblicz energię neutronu po śladzie.
        energy = neutron.energy - neutron.interaction_energy

        # Zaktualizuj stan reaktora.
        self.reactor.update_state(neutron, energy)

        return energy

    def calculate_neutron_scattering_cross_section(self, neutron, material):
        """
        Oblicza przekrój czynny rozpraszania neutronu na materiale.

        Argumenty:
            neutron (Neutron): Neutron do śledzenia.
            material (Material): Materiał, na którym neutron się rozprasza.

        Zwraca:
            float: Przekrój czynny rozpraszania neutronu na materiale.
        """

        # Oblicz energię neutronu.
        energy = neutron.energy

        # Oblicz gęstość neutronów w materiale.
        density = material.density

        # Oblicz masę atomową neutronu.
        atomic_mass = neutron.atomic_mass

        # Oblicz przekrój czynny rozpraszania neutronu na materiale.
        cross_section = density * atomic_mass * neutron.scattering_cross_section

        return cross_section

    def calculate_neutron_absorption_cross_section(self, neutron, material):
        """
        Oblicza przekrój czynny absorpcji neutronu na materiale.

        Argumenty:
            neutron (Neutron): Neutron do śledzenia.
            material (Material): Materiał, na którym neutron się absorbuje.

        Zwraca:
            float: Przekrój czynny absorpcji neutronu na materiale.
        """

        # Oblicz energię neutronu.
        energy = neutron.energy

        # Oblicz gęstość neutronów w materiale.
        density = material.density

        # Oblicz masę atomową neutronu.
        atomic_mass = neutron.atomic_mass

        # Oblicz przekrój czynny absorpcji neutronu na materiale.
        cross_section = density * atomic_mass * neutron.absorption_cross_section

        return cross_section


class ReactorController:

    def __init__(self, reactor):
        self.reactor = reactor

    def control_reactor(self, power_level):
        """
        Steruje reaktorem.

        Argumenty:
            power_level (float): Pożądany poziom mocy.

        Zwraca:
            bool: Czy udało się utrzymać pożądany poziom mocy?
        """

        # Oblicz aktualny poziom mocy.
        current_power_level = self.reactor.get_power_level()

        # Jeśli aktualny poziom mocy jest niższy od pożądanego poziomu mocy, zwiększ moc.
        if current_power_level < power_level:
            self.reactor.increase_power()

        # Jeśli aktualny poziom mocy jest wyższy od pożądanego poziomu mocy, zmniejsz moc.
        elif current_power_level > power_level:
            self.reactor.decrease_power()

        # Zwraca informację, czy udało się utrzymać pożądany poziom mocy.
        return current_power_level == power_level

    def get_current_power_level(self):
        """
        Pobiera aktualny poziom mocy reaktora.

        Zwraca:
            float: Aktualny poziom mocy reaktora.
        """

        return self.reactor.get_power_level()

    def get_maximum_power_level(self):
        """
        Pobiera maksymalny poziom mocy reaktora.

        Zwraca:
            float: Maksymalny poziom mocy reaktora.
        """

        return self.reactor.get_maximum_power_level()

    def get_minimum_power_level(self):
        """
        Pobiera minimalny poziom mocy reaktora.

        Zwraca:
            float: Minimalny poziom mocy reaktora.
        """

        return self.reactor.get_minimum_power_level()

    def increase_power(self, delta_power):
        """
        Zwiększa moc reaktora o daną wartość.

        Argumenty:
            delta_power (float): Wartość, o którą ma zostać zwiększona moc reaktora.
        """

        self.reactor.increase_power(delta_power)

    def decrease_power(self, delta_power):
        """
        Zmniejsza moc reaktora o daną wartość.

        Argumenty:
            delta_power (float): Wartość, o którą ma zostać zmniejszona moc reaktora.
        """

        self.reactor.decrease_power(delta_power)

    def get_status(self):
        """
        Pobiera status reaktora.

        Zwraca:
            str: Status reaktora.
        """

        return self.reactor.get_status()

    def get_errors(self):
        """
        Pobiera listę błędów reaktora.

        Zwraca:
            list: Lista błędów reaktora.
        """
        
        return self.reactor.get_errors()

    def clear_errors(self):
        """
        Czyści listę błędów reaktora.
        """

        self.reactor.clear_errors()


class RadiationDetector:

    def __init__(self, type, sensitivity):
        self.type = type
        self.sensitivity = sensitivity

    def detect_radiation(self, radiation):
        """
        Wykryj promieniowanie.

        Argumenty:
            radiation (Radiation): Promieniowanie do wykrycia.

        Zwraca:
            bool: Czy wykryto promieniowanie?
        """

        # Sprawdź, czy promieniowanie jest typu, którego detektor szuka.
        if self.type != radiation.type:
            return False

        # Sprawdź, czy promieniowanie jest na tyle silne, aby zostać wykryte.
        if radiation.intensity < self.sensitivity:
            return False

        # Detektor wykrył promieniowanie.
        return True

    def get_type(self):
        """
        Pobierz typ detektora promieniowania.

        Zwraca:
            str: Typ detektora promieniowania.
        """

        return self.type

    def get_sensitivity(self):
        """
        Pobierz czułość detektora promieniowania.

        Zwraca:
            float: Czułość detektora promieniowania.
        """

        return self.sensitivity

    def get_status(self):
        """
        Pobierz status detektora promieniowania.

        Zwraca:
            str: Status detektora promieniowania.
        """

        return self.status

    def get_errors(self):
        """
        Pobierz listę błędów detektora promieniowania.

        Zwraca:
            list: Lista błędów detektora promieniowania.
        """

        return self.errors

    def clear_errors(self):
        """
        Czyści listę błędów detektora promieniowania.
        """

        self.errors = []

    def set_sensitivity(self, sensitivity):
        """
        Ustaw czułość detektora promieniowania.

        Argumenty:
            sensitivity (float): Nowa czułość detektora promieniowania.
        """

        self.sensitivity = sensitivity

    def set_status(self, status):
        """
        Ustaw status detektora promieniowania.

        Argumenty:
            status (str): Nowy status detektora promieniowania.
        """

        self.status = status

    def add_error(self, error):
        """
        Dodaj błąd do listy błędów detektora promieniowania.

        Argumenty:
            error (str): Błąd do dodania.
        """

        self.errors.append(error)


class RadiationShielding:

    def __init__(self, material, thickness):
        self.material = material
        self.thickness = thickness

    def shield_radiation(self, radiation):
        """
        Chroń przed promieniowaniem.

        Argumenty:
            radiation (Radiation): Promieniowanie do ochrony.

        Zwraca:
            float: Redukcja dawki promieniowania.
        """

        # Oblicz redukcję dawki promieniowania.
        reduction = self.material.get_reduction_factor(radiation.type) * self.thickness

        # Zmniejsz dawkę promieniowania o wartość redukcji.
        radiation.dose -= reduction

        return reduction

    def get_material(self):
        """
        Pobierz materiał ekranu ochronnego.

        Zwraca:
            str: Materiał ekranu ochronnego.
        """

        return self.material

    def get_thickness(self):
        """
        Pobierz grubość ekranu ochronnego.

        Zwraca:
            float: Grubość ekranu ochronnego.
        """

        return self.thickness

    def get_status(self):
        """
        Pobierz status ekranu ochronnego.

        Zwraca:
            str: Status ekranu ochronnego.
        """

        return self.status

    def get_errors(self):
        """
        Pobierz listę błędów ekranu ochronnego.

        Zwraca:
            list: Lista błędów ekranu ochronnego.
        """

        return self.errors

    def clear_errors(self):
        """
        Czyści listę błędów ekranu ochronnego.
        """

        self.errors = []

    def set_status(self, status):
        """
        Ustaw status ekranu ochronnego.

        Argumenty:
            status (str): Nowy status ekranu ochronnego.
        """

        self.status = status

    def add_error(self, error):
        """
        Dodaj błąd do listy błędów ekranu ochronnego.

        Argumenty:
            error (str): Błąd do dodania.
        """

        self.errors.append(error)


class NuclearReactor:

    def __init__(self, type, fuel_material, shielding_material, power_level):
        self.type = type
        self.fuel_material = fuel_material
        self.shielding_material = shielding_material
        self.power_level = power_level

    def track_neutrons(self):
        """
        Śledzi neutrony w reaktorze.

        Zwraca:
            list: Lista neutronów w reaktorze.
        """

        # Utwórz listę neutronów.
        neutrons = []

        # Dodaj do listy neutrony z paliwa.
        for fuel_pellet in self.fuel:
            for neutron in fuel_pellet.neutrons:
                neutrons.append(neutron)

        # Dodaj do listy neutrony z osłony.
        for shielding_plate in self.shielding:
            for neutron in shielding_plate.neutrons:
                neutrons.append(neutron)

        return neutrons

    def control_reactor(self, power_level):
        """
        Kontroluje reaktor.

        Argumenty:
            power_level (float): Pożądany poziom mocy.

        Zwraca:
            bool: Czy udało się utrzymać pożądany poziom mocy?
        """

        # Oblicz aktualny poziom mocy.
        current_power_level = self.get_power_level()

        # Jeśli aktualny poziom mocy jest niższy od pożądanego poziomu mocy, zwiększ moc.
        if current_power_level < power_level:
            self.increase_power()

        # Jeśli aktualny poziom mocy jest wyższy od pożądanego poziomu mocy, zmniejsz moc.
        elif current_power_level > power_level:
            self.decrease_power()

        # Zwraca informację, czy udało się utrzymać pożądany poziom mocy.
        return current_power_level == power_level

    def detect_radiation(self, radiation):
        """
        Wykryj promieniowanie.

        Argumenty:
            radiation (Radiation): Promieniowanie do wykrycia.

        Zwraca:
            bool: Czy wykryto promieniowanie?
        """

        # Sprawdź, czy promieniowanie jest typu, którego detektor szuka.
        if self.detector.type != radiation.type:
            return False

        # Sprawdź, czy promieniowanie jest na tyle silne, aby zostać wykryte.
        if radiation.intensity < self.detector.sensitivity:
            return False

        # Detektor wykrył promieniowanie.
        return True

class NuclearPowerPlant:

    def __init__(self, type, number_of_reactors, power_output):
        self.type = type
        self.number_of_reactors = number_of_reactors
        self.power_output = power_output

    def generate_electricity(self):
        """
        Generuje energię elektryczną.

        Zwraca:
            float: Moc elektrowni.
        """

        # Oblicz moc elektrowni.
        power_output = self.number_of_reactors * self.reactor_power_output

        return power_output

    def dispose_of_nuclear_waste(self):
        """
        Usuwa odpady radioaktywne.

        Zwraca:
            bool: Czy udało się usunąć odpady radioaktywne?
        """

        # Sprawdź, czy jest dostępne składowisko odpadów radioaktywnych.
        if not self.waste_disposal_site:
            return False

        # Wyślij odpady radioaktywne do składowiska.
        self.waste_disposal_site.accept_waste(self.waste)

        return True

    def ensure_safety(self):
        """
        Zapewnia bezpieczeństwo elektrowni.

        Zwraca:
            bool: Czy udało się zapewnić bezpieczeństwo elektrowni?
        """

        # Sprawdź, czy reaktor jest w dobrym stanie technicznym.
        if not self.reactor.is_in_good_condition:
            return False

        # Sprawdź, czy osłona reaktora jest w dobrym stanie technicznym.
        if not self.shielding.is_in_good_condition:
            return False

        # Sprawdź, czy system bezpieczeństwa elektrowni działa poprawnie.
        if not self.safety_system.is_working_properly:
            return False

        return True

class NuclearWaste:

    def __init__(self, waste_type, half_life, toxicity, disposal_cost):
        self.waste_type = waste_type
        self.half_life = half_life
        self.toxicity = toxicity
        self.disposal_cost = disposal_cost

    def __str__(self):
        return f"{self.waste_type} ({self.half_life} years) ({self.toxicity}) ({self.disposal_cost})"

    def __repr__(self):
        return f"NuclearWaste({self.waste_type}, {self.half_life}, {self.toxicity}, {self.disposal_cost})"

    def get_waste_type(self):
        return self.waste_type

    def get_half_life(self):
        return self.half_life

    def get_toxicity(self):
        return self.toxicity

    def get_disposal_cost(self):
        return self.disposal_cost

       

class NuclearWasteDisposal:

    def __init__(self, location, waste_type, capacity):
        self.location = location
        self.waste_type = waste_type
        self.capacity = capacity

    def store_waste(self, waste):
        """
        Składuje odpady radioaktywne.

        Argumenty:
            waste (Waste): Odpady radioaktywne do składowania.

        Zwraca:
            bool: Czy udało się składować odpady radioaktywne?
        """

        # Sprawdź, czy składowisko jest dostępne.
        if not self.is_available:
            return False

        # Sprawdź, czy odpady są odpowiedniego typu.
        if self.waste_type != waste.type:
            return False

        # Sprawdź, czy składowisko jest pełne.
        if self.capacity <= 0:
            return False

        # Składuj odpady.
        self.waste.add_to_site(self)

        # Zmniejsz pojemność składowiska.
        self.capacity -= 1

        return True

    def monitor_site(self):
        """
        Monitoruje składowisko.

        Zwraca:
            dict: Stan składowiska.
        """

        # Zbierz informacje o stanie składowiska.
        status = {
            "location": self.location,
            "waste_type": self.waste_type,
            "capacity": self.capacity,
            "waste_count": len(self.waste),
        }

        return status

    def manage_site(self):
        """
        Zarządza składowiskiem.

        Zwraca:
            bool: Czy udało się zarządzać składowiskiem?
        """

        # Sprawdź, czy składowisko jest w dobrym stanie technicznym.
        if not self.is_in_good_condition:
            return False

        # Sprawdź, czy system monitorowania składowiska działa poprawnie.
        if not self.monitoring_system.is_working_properly:
            return False

        # Sprawdź, czy system bezpieczeństwa składowiska działa poprawnie.
        if not self.safety_system.is_working_properly:
            return False

        return True
    
    
    def calculate_activity(self):
        """
        Oblicza aktywność odpadu radioaktywnego.

        Zwraca:
            float: Aktywność odpadu radioaktywnego w becquerelach (Bq).
        """
        # Wzór na aktywność = lambda * N
        # gdzie lambda to stała rozpadu odpadu, a N to liczba rozpadów na sekundę.

        # Obliczenie stałej rozpadu odpadu (lambda) w 1/rok
        decay_constant = math.log(2) / self.half_life

        # Przeliczenie stałej rozpadu na 1/sekundę
        decay_constant_per_sec = decay_constant / (365 * 24 * 60 * 60)

        # Obliczenie liczby rozpadów na sekundę (aktywność)
        activity = decay_constant_per_sec * self.number_of_waste_items

        return activity

    def predict_disposal_cost(self, years):
        """
        Przewiduje koszt utylizacji odpadu radioaktywnego po określonym czasie.

        Argumenty:
            years (int): Liczba lat, po których przewidywany jest koszt utylizacji.

        Zwraca:
            float: Przewidywany koszt utylizacji odpadu radioaktywnego w jednostkach walutowych.
        """
        # Wzór na przewidywany koszt utylizacji = aktualny koszt * współczynnik inflacji^lata
        # Załóżmy, że koszt utylizacji rośnie o 3% rocznie.

        current_cost = self.disposal_cost
        inflation_rate = 0.03

        predicted_cost = current_cost * (1 + inflation_rate) ** years

        return predicted_cost

class NuclearSafety:

    def __init__(self, type, fuel_material, shielding_material, power_level):
        self.type = type
        self.fuel_material = fuel_material
        self.shielding_material = shielding_material
        self.power_level = power_level

    def detect_failure(self):
        """
        Wykryj awarię.

        Zwraca:
            bool: Czy wykryto awarię?
        """

        # Sprawdź, czy reaktor jest w stanie awarii.
        if self.reactor.is_in_failure_state:
            return True

        return False

    def prevent_failure(self):
        """
        Zapobiegaj awarii.

        Zwraca:
            bool: Czy udało się zapobiec awarii?
        """

        # Sprawdź, czy reaktor jest w stanie awarii.
        if self.reactor.is_in_failure_state:
            return False

        # Wyślij sygnał do systemu bezpieczeństwa reaktora.
        self.safety_system.alert()

        return True

    def respond_to_failure(self):
        """
        Odpowiedz na awarię.

        Zwraca:
            bool: Czy udało się odpowiedzieć na awarię?
        """

        # Sprawdź, czy reaktor jest w stanie awarii.
        if self.reactor.is_in_failure_state:
            return False

        # Wyłącz reaktor.
        self.reactor.shutdown()

        return True
class NuclearMaterial:

    def __init__(self, name, atomic_number, atomic_mass, half_life, fission_factor):
        self.name = name
        self.atomic_number = atomic_number
        self.atomic_mass = atomic_mass
        self.half_life = half_life
        self.fission_factor = fission_factor

    def __str__(self):
        return f"{self.name} ({self.atomic_number}, {self.atomic_mass})"

    def __repr__(self):
        return f"NuclearMaterial({self.name}, {self.atomic_number}, {self.atomic_mass}, {self.half_life}, {self.fission_factor})"

    def get_name(self):
        return self.name

    def get_atomic_number(self):
        return self.atomic_number

    def get_atomic_mass(self):
        return self.atomic_mass

    def get_half_life(self):
        return self.half_life

    def get_fission_factor(self):
        return self.fission_factor
    
class NuclearReaction:

    def __init__(self, reaction_type, energy_released, products):
        self.reaction_type = reaction_type
        self.energy_released = energy_released
        self.products = products

    def __str__(self):
        return f"{self.reaction_type} ({self.energy_released}) -> {self.products}"

    def __repr__(self):
        return f"NuclearReaction({self.reaction_type}, {self.energy_released}, {self.products})"

    def get_reaction_type(self):
        return self.reaction_type

    def get_energy_released(self):
        return self.energy_released

    def get_products(self):
        return self.products

class NuclearPower:

    def __init__(self, power_type, efficiency, cost, emissions):
        self.power_type = power_type
        self.efficiency = efficiency
        self.cost = cost
        self.emissions = emissions

    def __str__(self):
        return f"{self.power_type} ({self.efficiency}%) ({self.cost}) ({self.emissions})"

    def __repr__(self):
        return f"NuclearPower({self.power_type}, {self.efficiency}, {self.cost}, {self.emissions})"

    def get_power_type(self):
        return self.power_type

    def get_efficiency(self):
        return self.efficiency

    def get_cost(self):
        return self.cost

    def get_emissions(self):
        return self.emissions

    def calculate_power_output(self, input_energy):
        """
        Calculate the power output of the nuclear power plant.
        Input energy is the total energy supplied to the plant.
        """
        power_output = input_energy * self.efficiency / 100
        return power_output

    def compute_radiation_levels(self, fuel_mass):
        """
        Estimate the radiation levels based on the fuel mass and emissions.
        This is a simplified calculation and should not be used for real-world safety assessments.
        """
        radiation = fuel_mass * self.emissions
        return radiation

    def calculate_payback_period(self, annual_energy_output, initial_cost):
        """
        Calculate the payback period in years based on the annual energy output and initial cost.
        """
        if annual_energy_output <= 0:
            return None
        payback_period = initial_cost / (annual_energy_output * self.cost)
        return payback_period

    def calculate_roi(self, annual_energy_output, initial_cost, operation_years):
        """
        Calculate the return on investment (ROI) based on the annual energy output, initial cost, and operation years.
        """
        if annual_energy_output <= 0 or operation_years <= 0:
            return None
        total_revenue = annual_energy_output * self.cost * operation_years
        roi = (total_revenue - initial_cost) / initial_cost * 100
        return roi

    def calculate_lcoe(self, operation_years, annual_energy_output, discount_rate):
        """
        Calculate the levelized cost of electricity (LCOE) based on operation years, annual energy output, and discount rate.
        """
        if annual_energy_output <= 0 or operation_years <= 0 or discount_rate <= 0:
            return None
        present_value = self.cost * annual_energy_output / (1 - math.pow(1 + discount_rate, -operation_years))
        lcoe = present_value / (annual_energy_output * operation_years)
        return lcoe

class Atom:

    def __init__(self, number, mass, half_life, fission_factor):
        self.number = number
        self.mass = mass
        self.half_life = half_life
        self.fission_factor = fission_factor

    def __str__(self):
        return f"Atom: {self.number} {self.mass} {self.half_life} {self.fission_factor}"

    def __repr__(self):
        return f"Atom({self.number}, {self.mass}, {self.half_life}, {self.fission_factor})"

    def get_number(self):
        return self.number

    def get_mass(self):
        return self.mass

    def get_half_life(self):
        return self.half_life

    def get_fission_factor(self):
        return self.fission_factor

    def calculate_decay_rate(self):
        """
        Calculate the decay rate (disintegration rate) of the atom based on its half-life.
        """
        if self.half_life <= 0:
            return None
        decay_rate = math.log(2) / self.half_life
        return decay_rate

    def undergo_fission(self, neutron_count):
        """
        Check if the atom can undergo fission based on the number of neutrons supplied.
        If the neutron count is greater than the fission factor, fission occurs.
        """
        if neutron_count >= self.fission_factor:
            return True
        return False

    def calculate_fission_energy(self):
        """
        Calculate the energy released during fission based on the mass defect principle (E=mc^2).
        """
        c_squared = 299792458 ** 2  # Speed of light squared (m^2/s^2)
        mass_defect = self.mass - self.fission_factor * Atom.PROTON_MASS  # Assuming fission produces proton mass
        energy_released = mass_defect * c_squared
        return energy_released

class Molecule:

    def __init__(self, name, mass, chemical_formula, electron_configuration):
        self.name = name
        self.mass = mass
        self.chemical_formula = chemical_formula
        self.electron_configuration = electron_configuration

    def __str__(self):
        return f"Molecule: {self.name} {self.mass} {self.chemical_formula} {self.electron_configuration}"

    def __repr__(self):
        return f"Molecule({self.name}, {self.mass}, {self.chemical_formula}, {self.electron_configuration})"

    def get_name(self):
        return self.name

    def get_mass(self):
        return self.mass

    def get_chemical_formula(self):
        return self.chemical_formula

    def get_electron_configuration(self):
        return self.electron_configuration

    def count_atoms(self):
        """
        Count the number of atoms in the chemical formula of the molecule.
        Returns a dictionary with element symbols as keys and their corresponding counts as values.
        """
        atoms = {}
        current_element = ''
        count = ''
        for char in self.chemical_formula:
            if char.isalpha():
                if current_element:
                    if current_element in atoms:
                        atoms[current_element] += int(count) if count else 1
                    else:
                        atoms[current_element] = int(count) if count else 1
                current_element = char
                count = ''
            elif char.isdigit():
                count += char
        # Add the last element if it exists in the formula
        if current_element:
            if current_element in atoms:
                atoms[current_element] += int(count) if count else 1
            else:
                atoms[current_element] = int(count) if count else 1
        return atoms

    def calculate_molecular_weight(self):
        """
        Calculate the molecular weight of the molecule based on the atomic masses of its elements.
        Returns the total molecular weight in atomic mass units (u).
        """
        atoms = self.count_atoms()
        total_weight = 0.0
        for element, count in atoms.items():
            # Assuming you have a predefined dictionary 'atomic_masses' with atomic masses of elements.
            if element in atomic_masses:
                total_weight += atomic_masses[element] * count
        return total_weight

    def get_number_of_electrons(self):
        """
        Calculate the total number of electrons in the molecule based on its electron configuration.
        Returns the total number of electrons.
        """
        total_electrons = 0
        for orbital in self.electron_configuration:
            if orbital.isdigit():
                total_electrons += int(orbital)
        return total_electrons


class Crystal:

    def __init__(self, type, structure, crystallographic_parameters):
        self.type = type
        self.structure = structure
        self.crystallographic_parameters = crystallographic_parameters

    def __str__(self):
        return f"Crystal: {self.type} {self.structure} {self.crystallographic_parameters}"

    def __repr__(self):
        return f"Crystal({self.type}, {self.structure}, {self.crystallographic_parameters})"

    def get_type(self):
        return self.type

    def get_structure(self):
        return self.structure

    def get_crystallographic_parameters(self):
        return self.crystallographic_parameters

    def calculate_lattice_parameters(self):
        """
        Calculate the lattice parameters (a, b, c, alpha, beta, gamma) of the crystal.
        The crystallographic_parameters attribute should be a dictionary with keys 'a', 'b', 'c', 'alpha', 'beta', 'gamma'.
        """
        return self.crystallographic_parameters.get('a'), \
               self.crystallographic_parameters.get('b'), \
               self.crystallographic_parameters.get('c'), \
               self.crystallographic_parameters.get('alpha'), \
               self.crystallographic_parameters.get('beta'), \
               self.crystallographic_parameters.get('gamma')

    def calculate_density(self, mass):
        """
        Calculate the density of the crystal based on its mass and volume.
        """
        volume = self.calculate_volume()
        if volume is None or volume == 0:
            return None
        density = mass / volume
        return density

    def calculate_volume(self):
        """
        Calculate the volume of the crystal using the lattice parameters and crystallographic formulae.
        Returns the volume in cubic units.
        """
        a, b, c, alpha, beta, gamma = self.calculate_lattice_parameters()

        if None in [a, b, c, alpha, beta, gamma] or 0 in [a, b, c]:
            return None

        # Convert degrees to radians
        alpha_rad = math.radians(alpha)
        beta_rad = math.radians(beta)
        gamma_rad = math.radians(gamma)

        # Calculate volume using crystallographic formulae
        volume = a * b * c * math.sqrt(1 - math.cos(alpha_rad) ** 2 - math.cos(beta_rad) ** 2 - math.cos(gamma_rad) ** 2 + 2 * math.cos(alpha_rad) * math.cos(beta_rad) * math.cos(gamma_rad))

        return volume

class Solid:

    def __init__(self, melting_point, boiling_point, density, youngs_modulus):
        self.melting_point = melting_point
        self.boiling_point = boiling_point
        self.density = density
        self.youngs_modulus = youngs_modulus

    def __str__(self):
        return f"Solid: {self.melting_point} {self.boiling_point} {self.density} {self.youngs_modulus}"

    def __repr__(self):
        return f"Solid({self.melting_point}, {self.boiling_point}, {self.density}, {self.youngs_modulus})"

    def get_melting_point(self):
        return self.melting_point

    def get_boiling_point(self):
        return self.boiling_point

    def get_density(self):
        return self.density

    def get_youngs_modulus(self):
        return self.youngs_modulus

    def calculate_compressibility(self):
        """
        Calculate the compressibility of the solid material based on its Young's modulus.
        Compressibility is the reciprocal of the Young's modulus.
        """
        if self.youngs_modulus == 0:
            return None
        compressibility = 1 / self.youngs_modulus
        return compressibility

    def calculate_bulk_modulus(self):
        """
        Calculate the bulk modulus of the solid material based on its Young's modulus and density.
        """
        if self.youngs_modulus == 0 or self.density == 0:
            return None
        bulk_modulus = self.density * self.youngs_modulus
        return bulk_modulus

    def calculate_shear_modulus(self):
        """
        Calculate the shear modulus (modulus of rigidity) of the solid material based on its Young's modulus and Poisson's ratio.
        For simplicity, we assume Poisson's ratio to be 0.25.
        """
        if self.youngs_modulus == 0:
            return None
        poisson_ratio = 0.25
        shear_modulus = self.youngs_modulus / (2 * (1 + poisson_ratio))
        return shear_modulus

class Liquid:

    def __init__(self, melting_point, boiling_point, density, viscosity):
        self.melting_point = melting_point
        self.boiling_point = boiling_point
        self.density = density
        self.viscosity = viscosity

    def __str__(self):
        return f"Liquid: {self.melting_point} {self.boiling_point} {self.density} {self.viscosity}"

    def __repr__(self):
        return f"Liquid({self.melting_point}, {self.boiling_point}, {self.density}, {self.viscosity})"

    def get_melting_point(self):
        return self.melting_point

    def get_boiling_point(self):
        return self.boiling_point

    def get_density(self):
        return self.density

    def get_viscosity(self):
        return self.viscosity

    def calculate_dynamic_viscosity(self, flow_rate, area, length):
        """
        Calculate the dynamic viscosity of the liquid based on flow rate, area, and length.
        Dynamic viscosity (absolute viscosity) is measured in N·s/m² or Pa·s.
        """
        if flow_rate == 0 or area == 0 or length == 0:
            return None
        dynamic_viscosity = flow_rate * length / (area * 1000)  # Convert to Pa·s
        return dynamic_viscosity

    def calculate_kinematic_viscosity(self):
        """
        Calculate the kinematic viscosity of the liquid based on its dynamic viscosity and density.
        Kinematic viscosity is measured in m²/s.
        """
        if self.density == 0 or self.viscosity == 0:
            return None
        kinematic_viscosity = self.viscosity / self.density
        return kinematic_viscosity

    def calculate_reynolds_number(self, flow_rate, diameter):
        """
        Calculate the Reynolds number of the liquid flow based on flow rate and pipe diameter.
        Reynolds number is a dimensionless quantity used to predict flow regime (laminar or turbulent).
        """
        if flow_rate == 0 or diameter == 0 or self.viscosity == 0:
            return None
        reynolds_number = (flow_rate * diameter) / (self.viscosity * 1e-3)  # Convert viscosity to Pa·s
        return reynolds_number


class Gas:

    def __init__(self, melting_point, boiling_point, density, pressure):
        self.melting_point = melting_point
        self.boiling_point = boiling_point
        self.density = density
        self.pressure = pressure

    def __str__(self):
        return f"Gas: {self.melting_point} {self.boiling_point} {self.density} {self.pressure}"

    def __repr__(self):
        return f"Gas({self.melting_point}, {self.boiling_point}, {self.density}, {self.pressure})"

    def get_melting_point(self):
        return self.melting_point

    def get_boiling_point(self):
        return self.boiling_point

    def get_density(self):
        return self.density

    def get_pressure(self):
        return self.pressure

    def calculate_molar_volume(self, temperature, pressure):
        """
        Calculate the molar volume of the gas at a given temperature and pressure.
        Molar volume is the volume occupied by one mole of a gas at specific conditions.
        """
        if pressure <= 0 or temperature <= 0:
            return None

        # Ideal gas law constant (R) in L·atm/(mol·K)
        ideal_gas_constant = 0.0821

        # Convert pressure to atm and temperature to Kelvin
        pressure_atm = pressure / 101.325
        temperature_K = temperature + 273.15

        # Calculate molar volume using the ideal gas law: V = nRT/P
        molar_volume = (ideal_gas_constant * temperature_K) / pressure_atm
        return molar_volume

    def calculate_density_at_stp(self):
        """
        Calculate the density of the gas at standard temperature and pressure (STP).
        STP conditions are 0°C (273.15 K) and 1 atm (101.325 kPa).
        """
        if self.density <= 0:
            return None

        # Ideal gas law constant (R) in L·atm/(mol·K)
        ideal_gas_constant = 0.0821

        # Convert pressure to atm and temperature to Kelvin (STP conditions)
        pressure_atm = 1
        temperature_K = 273.15

        # Calculate the number of moles (n) of the gas using n = PV/RT
        moles = (self.density * 1000) / (ideal_gas_constant * temperature_K * pressure_atm)

        # Calculate molar volume using V = nRT/P
        molar_volume = (ideal_gas_constant * temperature_K) / pressure_atm

        # Calculate the gas density at STP using density = mass / volume
        gas_density_at_stp = moles / molar_volume
        return gas_density_at_stp

    def calculate_specific_gas_constant(self):
        """
        Calculate the specific gas constant (R_specific) for the gas.
        The specific gas constant is the universal gas constant (R) divided by the molar mass (M).
        """
        if self.density <= 0:
            return None

        # Universal gas constant (R) in L·atm/(mol·K)
        universal_gas_constant = 0.0821

        # Convert density to kg/m³ (density is in g/L)
        density_kg_m3 = self.density * 1000

        # Calculate molar volume at STP using the density at STP
        molar_volume_stp = self.calculate_molar_volume(0, 1)

        # Calculate the molar mass (M) using the ideal gas law: M = (density * V) / (RT)
        molar_mass = (density_kg_m3 * molar_volume_stp) / (universal_gas_constant * 273.15)

        # Calculate the specific gas constant (R_specific) using R_specific = R / M
        specific_gas_constant = universal_gas_constant / molar_mass
        return specific_gas_constant



class Wave:

    def __init__(self, wavelength, frequency, amplitude, velocity):
        self.wavelength = wavelength
        self.frequency = frequency
        self.amplitude = amplitude
        self.velocity = velocity

    def __str__(self):
        return f"Wave: {self.wavelength} {self.frequency} {self.amplitude} {self.velocity}"

    def __repr__(self):
        return f"Wave({self.wavelength}, {self.frequency}, {self.amplitude}, {self.velocity})"

    def get_wavelength(self):
        return self.wavelength

    def get_frequency(self):
        return self.frequency

    def get_amplitude(self):
        return self.amplitude

    def get_velocity(self):
        return self.velocity

    def calculate_wave_speed(self):
        """
        Calculate the wave speed (phase velocity) of the wave using the equation: v = λ * f
        """
        wave_speed = self.wavelength * self.frequency
        return wave_speed

    def calculate_period(self):
        """
        Calculate the period of the wave using the equation: T = 1 / f
        """
        if self.frequency == 0:
            return None
        period = 1 / self.frequency
        return period

    def calculate_energy(self):
        """
        Calculate the energy of the wave based on its amplitude and other properties.
        """
        # Assuming the wave is a transverse wave with linear displacement
        energy = 0.5 * (self.velocity ** 2) * (self.amplitude ** 2)
        return energy



class Particle:

    def __init__(self, mass, charge, velocity):
        self.mass = mass
        self.charge = charge
        self.velocity = velocity

    def __str__(self):
        return f"Particle: {self.mass} {self.charge} {self.velocity}"

    def __repr__(self):
        return f"Particle({self.mass}, {self.charge}, {self.velocity})"

    def get_mass(self):
        return self.mass

    def get_charge(self):
        return self.charge

    def get_velocity(self):
        return self.velocity

    def calculate_momentum(self):
        """
        Calculate the momentum of the particle using the formula: p = m * v
        """
        momentum = self.mass * self.velocity
        return momentum

    def calculate_kinetic_energy(self):
        """
        Calculate the kinetic energy of the particle using the formula: KE = 0.5 * m * v^2
        """
        kinetic_energy = 0.5 * self.mass * self.velocity ** 2
        return kinetic_energy

class Field:

    def __init__(self, intensity, direction, potential):
        self.intensity = intensity
        self.direction = direction
        self.potential = potential

    def __str__(self):
        return f"Field: {self.intensity} {self.direction} {self.potential}"

    def __repr__(self):
        return f"Field({self.intensity}, {self.direction}, {self.potential})"

    def get_intensity(self):
        return self.intensity

    def get_direction(self):
        return self.direction

    def get_potential(self):
        return self.potential

    def calculate_field_strength(self):
        """
        Calculate the field strength of the field using the intensity and direction.
        Field strength (E) is the magnitude of the intensity vector.
        """
        # Assuming intensity is a scalar value representing the magnitude of the vector
        field_strength = self.intensity
        return field_strength

    def calculate_energy_density(self):
        """
        Calculate the energy density of the field based on its intensity and potential.
        Energy density (u) of the field is given by u = 0.5 * ε * E^2,
        where ε is the permittivity of the medium and E is the field strength.
        For simplicity, we assume ε = 1 (vacuum).
        """
        field_strength = self.calculate_field_strength()
        if field_strength is None:
            return None

        # Permittivity of the vacuum (ε) is approximately 8.854 × 10^-12 F/m
        permittivity_vacuum = 8.854e-12

        energy_density = 0.5 * permittivity_vacuum * field_strength ** 2
        return energy_density



class Force:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Force: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Force({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_torque(self, lever_arm):
        """
        Calculate the torque (moment) of the force with respect to a given lever arm.
        Torque (τ) is given by τ = F * d, where F is the magnitude of the force and d is the lever arm.
        """
        if lever_arm <= 0:
            return None
        torque = self.magnitude * lever_arm
        return torque

    def calculate_work_done(self, displacement):
        """
        Calculate the work done by the force as it moves an object through a given displacement.
        Work done (W) is given by W = F * d * cos(θ),
        where F is the magnitude of the force, d is the displacement, and θ is the angle between the force and displacement vectors.
        """
        if displacement <= 0:
            return None

        # Calculate the angle between the force and displacement vectors (cosine of the angle)
        cos_theta = 1.0  # Assuming the force and displacement are aligned (θ = 0 degrees)

        work_done = self.magnitude * displacement * cos_theta
        return work_done


class Work:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Work: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Work({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_work_done(self, force_magnitude, displacement, angle=0):
        """
        Calculate the work done by a force as it moves an object through a given displacement.
        Work done (W) is given by W = F * d * cos(θ),
        where F is the magnitude of the force, d is the displacement, and θ is the angle between the force and displacement vectors.
        If the angle is not provided, it is assumed to be 0 degrees (force and displacement are aligned).
        """
        if displacement <= 0:
            return None

        # Calculate the angle between the force and displacement vectors (cosine of the angle)
        cos_theta = math.cos(math.radians(angle))

        work_done = force_magnitude * displacement * cos_theta
        return work_done

    def calculate_power(self, time):
        """
        Calculate the power associated with the work done over a given time.
        Power (P) is given by P = W / t, where W is the work done and t is the time taken.
        """
        if time <= 0:
            return None

        power = self.magnitude / time
        return power


class Energy:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Energy: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Energy({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_kinetic_energy(self, mass, velocity):
        """
        Calculate the kinetic energy of an object with a given mass and velocity.
        Kinetic energy (KE) is given by KE = 0.5 * m * v^2, where m is the mass and v is the velocity.
        """
        kinetic_energy = 0.5 * mass * velocity ** 2
        return kinetic_energy

    def calculate_potential_energy(self, mass, height, gravity=9.81):
        """
        Calculate the potential energy of an object with a given mass and height in a gravitational field.
        Potential energy (PE) is given by PE = m * g * h, where m is the mass, g is the acceleration due to gravity,
        and h is the height above a reference point.
        """
        potential_energy = mass * gravity * height
        return potential_energy

    def calculate_total_mechanical_energy(self, mass, height, velocity, gravity=9.81):
        """
        Calculate the total mechanical energy of an object in a gravitational field with a given mass, height, and velocity.
        Total mechanical energy (E) is the sum of kinetic energy and potential energy: E = KE + PE.
        """
        kinetic_energy = self.calculate_kinetic_energy(mass, velocity)
        potential_energy = self.calculate_potential_energy(mass, height, gravity)

        total_mechanical_energy = kinetic_energy + potential_energy
        return total_mechanical_energy



class Power:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Power: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Power({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_power(self, work_done, time):
        """
        Calculate the power given the work done and time taken.
        Power (P) is given by P = W / t, where W is the work done and t is the time taken.
        """
        if time <= 0:
            return None

        power = work_done / time
        return power



class Momentum:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Momentum: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Momentum({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_momentum(self, mass, velocity):
        """
        Calculate the momentum of an object with a given mass and velocity.
        Momentum (p) is given by p = m * v, where m is the mass and v is the velocity.
        """
        momentum = mass * velocity
        return momentum

class AngularMomentum:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Angular Momentum: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"AngularMomentum({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_angular_momentum(self, moment_of_inertia, angular_velocity):
        """
        Calculate the angular momentum of an object given its moment of inertia and angular velocity.
        Angular momentum (L) is given by L = I * ω, where I is the moment of inertia and ω is the angular velocity.
        """
        angular_momentum = moment_of_inertia * angular_velocity
        return angular_momentum



class Torque:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    
    def __str__(self):
        return f"Torque: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Torque({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_moment_of_inertia(self, mass, distance):
        return mass * distance ** 2

    def calculate_angular_acceleration(self, moment_of_inertia, torque):
        return torque / moment_of_inertia

    def calculate_angular_velocity(self, angular_acceleration, time):
        return angular_acceleration * time


class MomentOfInertia:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Moment Of Inertia: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Moment Of Inertia({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_angular_momentum(self, mass, velocity):
        return mass * velocity * self.direction

    def calculate_torque(self, angular_momentum):
        return angular_momentum / self.magnitude

    def calculate_rotational_kinetic_energy(self, angular_velocity):
        return 0.5 * self.magnitude * angular_velocity ** 2
    
class Entropy:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Entropy: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Entropy({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_information_entropy(self, probability):
        if probability == 0:
            return 0
        return -probability * math.log2(probability)

    def calculate_shannon_entropy(self, probabilities):
        if not probabilities:
            return 0
        return sum(-probability * math.log2(probability) for probability in probabilities)

def main():
    entropy = Entropy(10, math.pi / 2, (0, 0, 1))
    print(entropy)
    print(entropy.calculate_information_entropy(0.5))
    print(entropy.calculate_shannon_entropy([0.25, 0.25, 0.5]))
class Heat:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Heat: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Heat({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_thermal_energy(self, temperature_change):
        return self.magnitude * temperature_change

class Work:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Work: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Work({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_work(self, distance):
        return self.magnitude * distance

def main():
    heat = Heat(10, math.pi / 2, (0, 0, 1))
    print(heat)
    print(heat.calculate_thermal_energy(10))
    work = Work(10, math.pi / 2, (0, 0, 1))
    print(work)
    print(work.calculate_work(1))

if __name__ == "__main__":
    main()

class Work:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Work: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Work({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application


class Power:

    def __init__(self, magnitude, direction, point_of_application):
        self.magnitude = magnitude
        self.direction = direction
        self.point_of_application = point_of_application

    def __str__(self):
        return f"Power: {self.magnitude} {self.direction} {self.point_of_application}"

    def __repr__(self):
        return f"Power({self.magnitude}, {self.direction}, {self.point_of_application})"

    def get_magnitude(self):
        return self.magnitude

    def get_direction(self):
        return self.direction

    def get_point_of_application(self):
        return self.point_of_application

    def calculate_work(self, distance):
        return self.magnitude * distance

    def calculate_force(self, distance):
        return self.magnitude / distance

def main():
    power = Power(10, math.pi / 2, (0, 0, 1))
    print(power)
    print(power.calculate_work(1))
    print(power.calculate_force(1))
    
class Atom:

    def __init__(self, name, atomic_number, atomic_mass, electron_configuration, ionization_energy, electronegativity):
        self.name = name
        self.atomic_number = atomic_number
        self.atomic_mass = atomic_mass
        self.electron_configuration = electron_configuration
        self.ionization_energy = ionization_energy
        self.electronegativity = electronegativity

    def __str__(self):
        return f"Atom: {self.name} {self.atomic_number} {self.atomic_mass} {self.electron_configuration} {self.ionization_energy} {self.electronegativity}"

    def __repr__(self):
        return f"Atom({self.name}, {self.atomic_number}, {self.atomic_mass}, {self.electron_configuration}, {self.ionization_energy}, {self.electronegativity})"

    def calculate_ionization_potential(self):
        return math.sqrt(self.ionization_energy)

    def calculate_electron_affinity(self):
        return -math.sqrt(self.electronegativity)

    def print_atom(self):
        print(f"Atom: {self.name} {self.atomic_number} {self.atomic_mass}")
        print(f"  Electron configuration: {self.electron_configuration}")
        print(f"  Ionization energy: {self.ionization_energy}")
        print(f"  Electronegativity: {self.electronegativity}")

def main():
    atom = Atom("Hydrogen", 1, 1.00794, "1s1", 13.602, 2.20)
    print(atom)
    print(atom.calculate_ionization_potential())
    print(atom.calculate_electron_affinity())
    atom.print_atom()

class Bond:

    def __init__(self, bond_type, bond_order, bond_length, bond_angle, bond_dipole_moment):
        self.bond_type = bond_type
        self.bond_order = bond_order
        self.bond_length = bond_length
        self.bond_angle = bond_angle
        self.bond_dipole_moment = bond_dipole_moment

    def __str__(self):
        return f"Bond: {self.bond_type} {self.bond_order} {self.bond_length} {self.bond_angle} {self.bond_dipole_moment}"

    def __repr__(self):
        return f"Bond({self.bond_type}, {self.bond_order}, {self.bond_length}, {self.bond_angle}, {self.bond_dipole_moment})"

    def calculate_bond_energy(self):
        if self.bond_type == "single":
            return 413.4
        elif self.bond_type == "double":
            return 611.0
        elif self.bond_type == "triple":
            return 836.8
        else:
            raise ValueError(f"Unknown bond type: {self.bond_type}")

    def calculate_force_constant(self):
        if self.bond_type == "single":
            return 510.0
        elif self.bond_type == "double":
            return 2620.0
        elif self.bond_type == "triple":
            return 8360.0
        else:
            raise ValueError(f"Unknown bond type: {self.bond_type}")
class Molecule:

    def __init__(self, name, formula, atoms, bonds):
        self.name = name
        self.formula = formula
        self.atoms = atoms
        self.bonds = bonds

    def __str__(self):
        return f"Molecule: {self.name} {self.formula} {self.atoms} {self.bonds}"

    def __repr__(self):
        return f"Molecule({self.name}, {self.formula}, {self.atoms}, {self.bonds})"

    def calculate_molecular_weight(self):
        return sum(atom.atomic_mass for atom in self.atoms)

    def calculate_bond_angles(self):
        return [angle for angle in self.bonds]

    def calculate_dipole_moment(self):
        return sum(atom.charge * atom.position for atom in self.atoms)

    def print_molecule(self):
        print(f"Molecule: {self.name} {self.formula}")
        for atom in self.atoms:
            print(f"  {atom}")
        for bond in self.bonds:
            print(f"  {bond}")

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

class Weapon:

    def __init__(self, type, power, range):
        self.type = type
        self.power = power
        self.range = range

    def __str__(self):
        return f"Weapon type: {self.type} | Power: {self.power} | Range: {self.range}"

    def __repr__(self):
        return f"Weapon({self.type}, {self.power}, {self.range})"

    def get_type(self):
        return self.type

    def get_power(self):
        return self.power

    def get_range(self):
        return self.range

    def set_type(self, type):
        self.type = type

    def set_power(self, power):
        self.power = power

    def set_range(self, range):
        self.range = range

    def predict_damage(self, target):
        """Predicts the damage that the weapon will do to the target.

        Args:
            target (Weapon): The target of the attack.

        Returns:
            float: The predicted damage.
        """

        # Load the training data.
        data = pd.read_csv("weapons.csv")

        # Create the features and target variables.
        features = data[["type", "power", "range"]]
        target = data["damage"]

        # Create the linear regression model.
        model = LinearRegression()

        # Fit the model to the training data.
        model.fit(features, target)

        # Predict the damage for the target weapon.
        prediction = model.predict(features.loc[features["type"] == self.type, :])

        return prediction

    def predict_probability_of_hit(self, target):
        """Predicts the probability that the weapon will hit the target.

        Args:
            target (Weapon): The target of the attack.

        Returns:
            float: The predicted probability of a hit.
        """

        # Load the training data.
        data = pd.read_csv("weapons.csv")

        # Create the features and target variables.
        features = data[["type", "power", "range"]]
        target = data["hit"]

        # Create the logistic regression model.
        model = LogisticRegression()

        # Fit the model to the training data.
        model.fit(features, target)

        # Predict the probability of a hit for the target weapon.
        prediction = model.predict_proba(features.loc[features["type"] == self.type, :])

        return prediction[0][1]

    def predict_damage_and_probability_of_hit(self, target):
        """Predicts the damage that the weapon will do to the target and the probability of a hit.

        Args:
            target (Weapon): The target of the attack.

        Returns:
            tuple: The predicted damage and probability of a hit.
        """

        damage = self.predict_damage(target)
        probability_of_hit = self.predict_probability_of_hit(target)

        return damage, probability_of_hit

    def get_most_likely_target(self, targets):
        """Gets the most likely target for the weapon.

        Args:
            targets (list): A list of targets.

        Returns:
            Weapon: The most likely target.
        """

        # Calculate the probability of hitting each target.
        probabilities = [self.predict_probability_of_hit(target) for target in targets]

        # Get the target with the highest probability.
        most_likely_target = targets[np.argmax(probabilities)]

        # If the probability of hitting the most likely target is less than 0.5,
        # return None.
        if probabilities[np.argmax(probabilities)] < 0.5:
            return None

        return most_likely_target

    def get_expected_damage(self, targets):
        """Gets the expected damage that the weapon will do to the targets.

        Args:
            targets (list): A list of targets.

        Returns:
            float: The expected damage.
        """

        # Calculate the expected damage for each target.
        expected_damage = []
        for target in targets:
            damage = self.predict_damage(target)
            probability = self.predict_probability_of_hit(target)
            expected_damage.append(damage * probability)

        # Return the sum of the expected damage.
        return sum(expected_damage)

if __name__ == "__main__":
    fission_analysis = FissionAnalysis()
    start_mass = 235
    end_mass = 250
    masses, energies = fission_analysis.collect_data(start_mass, end_mass)
    fission_analysis.plot_energy_vs_mass(masses, energies)
    torque = Torque(10, math.pi / 2, (0, 0, 1))
    print(torque)
    print(torque.calculate_moment_of_inertia(1, 1))
    print(torque.calculate_angular_acceleration(1, 10))
    print(torque.calculate_angular_velocity(math.pi / 2, 1))
    moment_of_inertia = MomentOfInertia(1, math.pi / 2, (0, 0, 1))
    molecule = Molecule("Water", "H2O", ["H", "H", "O"], ["O-H", "H-O"])
    print(molecule)
    print(molecule.calculate_molecular_weight())
    print(molecule.calculate_bond_angles())
    print(molecule.calculate_dipole_moment())
    molecule.print_molecule()
    print(moment_of_inertia)
    print(moment_of_inertia.calculate_angular_momentum(1, 1))
    print(moment_of_inertia.calculate_torque(1))
    print(moment_of_inertia.calculate_rotational_kinetic_energy(math.pi / 2))
    entropy = Entropy(1, math.pi / 2, (0, 0, 1))
    print(entropy)
    print(entropy.calculate_information_entropy(0.5))
    print(entropy.calculate_shannon_entropy([0.5, 0.5]))
    bond = Bond("single", 1, 0.12, 104.5, 1.6)
    print(bond)
    print(bond.calculate_bond_energy())
    print(bond.calculate_force_constant())
    weapon = Weapon("gun", 10, 100)

    # Print the weapon's information.
    print(weapon)

    # Predict the damage that the weapon will do to a target.
    target = Weapon("tank", 1000, 10000)
    damage = weapon.predict_damage(target)
    print("The weapon will do {} damage to the target.".format(damage))

    # Predict the probability that the weapon will hit a target.
    probability = weapon.predict_probability_of_hit(target)
    print("The weapon has a {}% chance of hitting the target.".format(probability))

    # Get the most likely target for the weapon.
    targets = [Weapon("tank", 1000, 10000), Weapon("plane", 10000, 100000)]
    most_likely_target = weapon.get_most_likely_target(targets)
    print("The most likely target for the weapon is {}.".format(most_likely_target))

    # Get the expected damage that the weapon will do to the targets.
    expected_damage = weapon.get_expected_damage(targets)
    print(f"The weapon will do {expected_damage} expected damage to the targets.")


class WeaponTest(unittest.TestCase):

    def test_predict_damage(self):
        weapon = Weapon("gun", 10, 100)
        target = Weapon("tank", 1000, 10000)
        damage = weapon.predict_damage(target)
        self.assertEqual(damage, 1000)

    def test_predict_probability_of_hit(self):
        weapon = Weapon("gun", 10, 100)
        target = Weapon("tank", 1000, 10000)
        probability = weapon.predict_probability_of_hit(target)
        self.assertEqual(probability, 0.5)

    def test_get_most_likely_target(self):
        weapon = Weapon("gun", 10, 100)
        targets = [Weapon("tank", 1000, 10000), Weapon("plane", 10000, 100000)]
        most_likely_target = weapon.get_most_likely_target(targets)
        self.assertEqual(most_likely_target, targets[0])

    def test_get_expected_damage(self):
        weapon = Weapon("gun", 10, 100)
        targets = [Weapon("tank", 1000, 10000), Weapon("plane", 10000, 100000)]
        expected_damage = weapon.get_expected_damage(targets)
        self.assertEqual(expected_damage, 5000)
    
    def test_predict_damage_with_armor(self):
        weapon = Weapon("gun", 10, 100)
        target = Weapon("tank", 1000, 10000, armor=100)
        damage = weapon.predict_damage(target)
        self.assertEqual(damage, 900)

    def test_predict_probability_of_hit_with_health(self):
        weapon = Weapon("gun", 10, 100)
        target = Weapon("tank", 1000, 10000, health=500)
        probability = weapon.predict_probability_of_hit(target)
        self.assertEqual(probability, 0.75)

    def test_get_most_likely_target_with_multiple_targets(self):
        weapon = Weapon("gun", 10, 100)
        targets = [
            Weapon("tank", 1000, 10000),
            Weapon("plane", 10000, 100000),
            Weapon("ship", 100000, 1000000),
        ]
        most_likely_target = weapon.get_most_likely_target(targets)
        self.assertEqual(most_likely_target, targets[0])

    def test_get_expected_damage_with_multiple_targets(self):
        weapon = Weapon("gun", 10, 100)
        targets = [
            Weapon("tank", 1000, 10000),
            Weapon("plane", 10000, 100000),
            Weapon("ship", 100000, 1000000),
        ]
        expected_damage = weapon.get_expected_damage(targets)
        self.assertEqual(expected_damage, 35000)
        
    
class PhysicsTest(unittest.TestCase):

    def test_calculate_power(self):
        power = Power(100, "forward", "center")
        work_done = 1000
        time = 10
        calculated_power = power.calculate_power(work_done, time)
        self.assertEqual(calculated_power, 100)

    def test_calculate_momentum(self):
        momentum = Momentum(100, "forward", "center")
        mass = 10
        velocity = 10
        calculated_momentum = momentum.calculate_momentum(mass, velocity)
        self.assertEqual(calculated_momentum, 1000)

    def test_calculate_angular_momentum(self):
        angular_momentum = AngularMomentum(100, "forward", "center")
        moment_of_inertia = 10
        angular_velocity = 10
        calculated_angular_momentum = angular_momentum.calculate_angular_momentum(moment_of_inertia, angular_velocity)
        self.assertEqual(calculated_angular_momentum, 1000)

    def test_calculate_torque(self):
        torque = Torque(100, "forward", "center")
        moment_of_inertia = 10
        angular_momentum = 100
        calculated_torque = torque.calculate_torque(angular_momentum)
        self.assertEqual(calculated_torque, 10)

    def test_calculate_moment_of_inertia(self):
        moment_of_inertia = MomentOfInertia(100, "forward", "center")
        mass = 10
        distance = 10
        calculated_moment_of_inertia = moment_of_inertia.calculate_moment_of_inertia(mass, distance)
        self.assertEqual(calculated_moment_of_inertia, 1000)


if __name__ == "__main__":
    unittest.main()
import unittest

from nuclear import NuclearWasteDisposal, NuclearSafety, NuclearMaterial, NuclearReaction, NuclearPower

class NuclearWasteDisposalTest(unittest.TestCase):

    def test_store_waste(self):
        disposal = NuclearWasteDisposal("Location", "Type", 10)
        waste = NuclearMaterial("Waste", 1, 1, 1, 1)

        # Test storing waste when the disposal is available.
        self.assertTrue(disposal.store_waste(waste))

        # Test storing waste when the disposal is not available.
        disposal.is_available = False
        self.assertFalse(disposal.store_waste(waste))

        # Test storing waste when the waste is not of the correct type.
        waste = NuclearMaterial("Waste", 1, 1, 1, 2)
        self.assertFalse(disposal.store_waste(waste))

        # Test storing waste when the disposal is full.
        for _ in range(10):
            disposal.store_waste(waste)
        self.assertFalse(disposal.store_waste(waste))

    def test_monitor_site(self):
        disposal = NuclearWasteDisposal("Location", "Type", 10)

        # Test monitoring the site when the disposal is available.
        status = disposal.monitor_site()
        self.assertEqual(status["location"], "Location")
        self.assertEqual(status["waste_type"], "Type")
        self.assertEqual(status["capacity"], 10)
        self.assertEqual(status["waste_count"], 0)

        # Test monitoring the site when the disposal is not available.
        disposal.is_available = False
        status = disposal.monitor_site()
        self.assertIsNone(status)

    def test_manage_site(self):
        disposal = NuclearWasteDisposal("Location", "Type", 10)

        # Test managing the site when the disposal is available.
        self.assertTrue(disposal.manage_site())

        # Test managing the site when the disposal is not available.
        disposal.is_available = False
        self.assertFalse(disposal.manage_site())

class NuclearSafetyTest(unittest.TestCase):

    def test_detect_failure(self):
        safety = NuclearSafety("Type", "Fuel", "Shielding", 100)

        # Test detecting a failure when the reactor is in a failure state.
        safety.reactor.is_in_failure_state = True
        self.assertTrue(safety.detect_failure())

        # Test detecting a failure when the reactor is not in a failure state.
        safety.reactor.is_in_failure_state = False
        self.assertFalse(safety.detect_failure())

    def test_prevent_failure(self):
        safety = NuclearSafety("Type", "Fuel", "Shielding", 100)

        # Test preventing a failure when the reactor is in a failure state.
        safety.reactor.is_in_failure_state = True
        self.assertTrue(safety.prevent_failure())

        # Test preventing a failure when the reactor is not in a failure state.
        safety.reactor.is_in_failure_state = False
        self.assertFalse(safety.prevent_failure())

    def test_respond_to_failure(self):
        safety = NuclearSafety("Type", "Fuel", "Shielding", 100)

        # Test responding to a failure when the reactor is in a failure state.
        safety.reactor.is_in_failure_state = True
        self.assertTrue(safety.respond_to_failure())

        # Test responding to a failure when the reactor is not in a failure state.
        safety.reactor.is_in_failure_state = False
        self.assertFalse(safety.respond_to_failure())

class NuclearMaterialTest(unittest.TestCase):

    def test_init(self):
        material = NuclearMaterial("Material", 1, 1, 1, 1)

        self.assertEqual(material.name, "Material")
        self.assertEqual(material.atomic_number, 1)
        self.assertEqual(material.atomic_mass, 1)
        self.assertEqual(material.half_life, 1)
        self.assertEqual(material.fission_factor, 1)
