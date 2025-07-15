#ifndef BMP_FDCT_c
#define BMP_FDCT_c

#include <tuple>
#include <string>
#include <iostream>
#include <filesystem>
#include <cstdlib> // per std::stoi
#include <opencv2/imgcodecs.hpp>
#include <opencv2/opencv.hpp>
#include "tinyfiledialogs.h"

class BMPFDCT
{
public:
    BMPFDCT() = default;
    ~BMPFDCT() = default;

    void test_it()
    {
        std::tuple<std::string, int, int> result = ask_user_for_input();
        std::string bmp = std::get<0>(result);
        int f = std::get<1>(result);
        int d = std::get<2>(result);
        cv::Mat image = cv::imread(bmp, cv::IMREAD_GRAYSCALE);
        if (image.empty()) {
            std::cerr << "Errore nel caricamento dell'immagine\n";
            return;
        }

        process_and_display_dct(image, f, d, bmp);
    }

    void process_and_display_dct(const cv::Mat &image, int f, int d, const std::string& bmp_path)
    {
        // Controllo dimensioni immagine e scarto gli avanzi (bordo)
        int width = image.cols - (image.cols % f);
        int height = image.rows - (image.rows % f);
        // Matrice risultato, inizializzata a zero, dimensione adattata
        cv::Mat result_image = cv::Mat::zeros(height, width, CV_8UC1);

        // Ciclo sui blocchi quadrati f x f
        for (int row = 0; row < height; row += f)
        {
            for (int col = 0; col < width; col += f)
            {
                // Estraggo il blocco di dimensione f x f dall'immagine originale
                cv::Mat block = image(cv::Rect(col, row, f, f));

                // Converto il blocco in double (CV_64F) per la DCT
                cv::Mat block_double;
                block.convertTo(block_double, CV_64F);

                // Calcolo la DCT2 del blocco
                cv::Mat dct_block;
                cv::dct(block_double, dct_block);

                // Applico la soglia: azzero i coefficienti dove k + l >= d
                // k = riga, l = colonna nella matrice dei coefficienti DCT
                for (int k = 0; k < f; ++k)
                {
                    for (int l = 0; l < f; ++l)
                    {
                        if ((k + l) >= d)
                        {
                            dct_block.at<double>(k, l) = 0.0;
                        }
                    }
                }

                // Calcolo l'IDCT2 per ottenere il blocco ricostruito
                cv::Mat idct_block;
                cv::idct(dct_block, idct_block);

                // Arrotondo i valori e li limito tra 0 e 255 per formare un'immagine valida
                for (int i = 0; i < f; ++i)
                {
                    for (int j = 0; j < f; ++j)
                    {
                        double val = std::round(idct_block.at<double>(i, j));
                        if (val < 0)
                            val = 0;
                        else if (val > 255)
                            val = 255;
                        idct_block.at<double>(i, j) = val;
                    }
                }

                // Converto il blocco da double a uchar (CV_8U)
                cv::Mat block_uchar;
                idct_block.convertTo(block_uchar, CV_8U);

                // Copio il blocco ricostruito nella posizione corretta dell'immagine risultato
                block_uchar.copyTo(result_image(cv::Rect(col, row, f, f)));
            }
        }

        // Visualizzo affiancate l'immagine originale e quella ricostruita
        cv::Mat display;
        // Crop l’immagine originale per farla della stessa dimensione del risultato (scarto avanzi)
        cv::Mat original_cropped = image(cv::Rect(0, 0, width, height));

        //int height = original_cropped.rows;
        int sep_width = 5;                    // spessore del separatore in pixel
        cv::Mat separator = cv::Mat::zeros(height, sep_width, original_cropped.type());

        // Concatena originale | separatore | DCT
        std::vector<cv::Mat> parts = { original_cropped, separator, result_image };
        cv::hconcat(parts, display);

        cv::imshow("Originale (sinistra) e Modificata (destra)", display);
        cv::waitKey(0);

        std::string save_path = "immagini/after_dct/";
        //Costruisci il nome del file
        std::string base_name = std::filesystem::path(bmp_path).stem().string();  // es: "20x20"
        std::string output_file = save_path + base_name + " F= " + std::to_string(f) + " d= " + std::to_string(d) + "_dct.bmp";

        //Salva immagine modificata
        cv::imwrite(output_file, result_image);
    }

private:
    std::tuple<std::string, int, int> ask_user_for_input()
    {
        const char *filters[] = {"*.bmp"};

        const char *bmp_path = tinyfd_openFileDialog(
            "Seleziona un file BMP in toni di grigio",
            "",
            1,
            filters,
            "Immagini BMP",
            0);

        if (!bmp_path)
        {
            std::cerr << "Nessun file selezionato. Uscita.\n";
            return std::make_tuple("", -1, -1);
        }

        const char *f_str = tinyfd_inputBox(
            "Inserisci F",
            "Ampiezza delle finestrelle (intero positivo):",
            "8");

        if (!f_str)
        {
            std::cerr << "Input F annullato. Uscita.\n";
            return std::make_tuple("", -1, -1);
        }

        int F = std::stoi(f_str);

        const char *d_str = tinyfd_inputBox(
            "Inserisci d",
            "Soglia di taglio delle frequenze (intero tra 0 e 2F-2):",
            "5");

        if (!d_str)
        {
            std::cerr << "Input d annullato. Uscita.\n";
            return std::make_tuple("", -1, -1);
        }

        int d = std::stoi(d_str);

        if (F <= 0)
        {
            std::cerr << "F deve essere positivo.\n";
            return std::make_tuple("", -1, -1);
        }

        if (d < 0 || d > 2 * F - 2)
        {
            std::cerr << "d deve essere compreso tra 0 e " << 2 * F - 2 << ".\n";
            return std::make_tuple("", -1, -1);
        }
        std::cout << std::string(bmp_path) << " " << F << " " << d << std::endl;
        return std::make_tuple(std::string(bmp_path), F, d);
    }
};

#endif
