// blueseltzer (c) 2025 orthopteroid@gmail.com
// MIT license

///////////////////
// platform headers

#if defined(WIN32)

#include <windows.h>
#include <comdef.h>
#include <combaseapi.h>
#include <wincodec.h>
#include <shlwapi.h>

#elif defined(LINUX)

#include <stdio.h>
#include <memory.h>
#include <setjmp.h>
#include <jpeglib.h>
#include <jerror.h>

#endif // platform headers

///////////////////
// common headers and basic types

#include <memory>
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <numbers>

using namespace std;

struct surf_t
{
    uint32_t width, height;
    uint8_t* data;
};

struct point_t {
    int x, y;
};

struct bbox_t {
    int xmin, xmax, ymin, ymax;
};

struct Surface;

///////////////////
// platform dependent stuff

#if defined(WIN32)

struct Image
{
    const wchar_t* wzFilename;
    HRESULT hr;
    IWICImagingFactory* pIWICFactory;
    IWICBitmapDecoder* pIDecoder;
    IWICBitmapFrameDecode* pIFrame;
    IWICFormatConverter* pIConverter;

    Image()
    {
        wzFilename = NULL;
        hr = S_OK;
        pIWICFactory = NULL;
        pIDecoder = NULL;
        pIFrame = NULL;
        pIConverter = NULL;
    }
    virtual ~Image()
    {
        if (pIFrame) pIFrame->Release();
        if (pIConverter) pIConverter->Release();
        if (pIDecoder) pIDecoder->Release();
        if (pIWICFactory) pIWICFactory->Release();
        CoUninitialize();
    }
    
    void SetFile(const wchar_t* wz)
    {
        wzFilename = wz;
    }

    // https://gist.github.com/xebecnan/6d070c93fb69f40c3673
    std::string GetFile()
    {
        if (!wzFilename) return std::string();
        size_t wclen = wcslen(wzFilename);
        size_t mblen = WideCharToMultiByte(CP_UTF8, 0, wzFilename, wclen, 0, 0, NULL, NULL);
        char* buff = new char[mblen + 1];
        WideCharToMultiByte(CP_UTF8, 0, wzFilename, wclen, buff, mblen, NULL, NULL);
        buff[mblen] = '\0';
        std::string str(buff);
        delete []buff;
        return str;
    }

    void Load()
    {
        const char* szContext = NULL;

        szContext = "CoInitializeEx";
        hr = CoInitializeEx(NULL, COINIT_MULTITHREADED);
        if (FAILED(hr)) goto err;

        szContext = "CoCreateInstance";
        hr = CoCreateInstance(
            CLSID_WICImagingFactory,
            NULL,
            CLSCTX_INPROC_SERVER,
            IID_PPV_ARGS(&pIWICFactory)
        );
        if (FAILED(hr)) goto err;

        szContext = "IWICImagingFactory::CreateDecoderFromFilename";
        hr = pIWICFactory->CreateDecoderFromFilename(
            wzFilename,
            NULL,                           // Do not prefer a particular vendor
            GENERIC_READ,                   // Desired read access to the file
            WICDecodeMetadataCacheOnDemand, // Cache metadata when needed
            &pIDecoder                      // Pointer to the decoder
        );
        if (FAILED(hr)) goto err;

        szContext = "IWICImagingFactory::CreateFormatConverter";
        hr = pIWICFactory->CreateFormatConverter(&pIConverter);
        if (FAILED(hr)) goto err;

        // rescale?
        // https://github.com/MicrosoftDocs/win32/blob/docs/desktop-src/wic/-wic-bitmapsources-howto-scale.md

        szContext = "IWICBitmapDecoder::GetFrame";
        hr = pIDecoder->GetFrame(0, &pIFrame);
        if (FAILED(hr)) goto err;

        szContext = "IWICFormatConverter::Initialize";
        hr = pIConverter->Initialize(
            pIFrame,
            GUID_WICPixelFormat8bppIndexed,
            WICBitmapDitherTypeNone,
            0, 0,
            WICBitmapPaletteTypeFixedGray256
        );
        if (FAILED(hr)) goto err;

        return;
    err:
        _com_error err(hr);
        cout << "(win32) " << szContext << ": " << err.ErrorMessage() << endl;
        exit(1);
    }

    unique_ptr<Surface> GetSurface()
    {
        const char* szContext = NULL;
        surf_t surf;

        szContext = "CreateFormatConverter::GetSize";
        hr = pIConverter->GetSize(&surf.width, &surf.height);
        if (FAILED(hr)) goto err;

        surf.data = new uint8_t[surf.width * surf.height](0);

        szContext = "CreateFormatConverter::CopyPixels";
        hr = pIConverter->CopyPixels(NULL, surf.width, surf.width * surf.height, surf.data);
        if (FAILED(hr)) goto err;

        return make_unique<Surface>(surf);
    err:
        _com_error err(hr);
        cout << "(win32) " << szContext << ": " << err.ErrorMessage() << endl;
        exit(1);
    }
};

void app(Image& image); // fwd declaration

// win32 console app
int wmain(int argc, wchar_t* argv[])
{
    Image image;
    image.SetFile(argv[1]);
    app(image);
    return 0;
}

#elif defined(LINUX)

// https://www.tspi.at/2020/03/20/libjpegexample.html#gsc.tab=0
struct Image
{
    char* szFilename;
    surf_t surf;

    Image()
    {
    }
    virtual ~Image()
    {
    }

    void SetFile(char* sz)
    {
        szFilename = sz;
    }

    std::string GetFile()
    {
        return std::string(szFilename);
    }

    void Load()
    {
        const char* szContext = "";
        FILE* fHandle = NULL;
        uint8_t* pBuffer = 0;

        struct jpeg_decompress_struct info;
        struct jpeg_error_mgr err;

        szContext = "fopen";
        if ((fHandle = fopen(szFilename, "rb")) == NULL) goto err;

        info.err = jpeg_std_error(&err);
        jpeg_create_decompress(&info);
        jpeg_stdio_src(&info, fHandle);

        szContext = "jpeg_read_header";
        if(jpeg_read_header(&info, TRUE) == 0) goto err;

        szContext = "jpeg_start_decompress";
        if (jpeg_start_decompress(&info) == false) goto err;

        szContext = "info.num_components != 3";
        if(info.num_components != 3) goto err;

        surf.width = info.output_width;
        surf.height=info.output_height;
        surf.data = new uint8_t[info.output_width * info.output_height];
        pBuffer = new uint8_t[info.output_width * info.num_components];

        while (info.output_scanline < info.output_height)
        {
            uint8_t* ppBuffer[1];
            ppBuffer[0] = pBuffer;
            szContext = "jpeg_read_scanlines";
            JDIMENSION n = jpeg_read_scanlines(&info, ppBuffer, 1);
            if (n!=1) goto err;
            for(int i=0; i<info.output_width; i++)
                surf.data[info.output_width * (info.output_scanline-1) + i] =
                    (uint8_t)(0.299f * (float)pBuffer[i*3+0] + 0.587f * (float)pBuffer[i*3+1] + 0.114f * (float)pBuffer[i*3+2]);
        }

        szContext = "jpeg_finish_decompress";
        if (jpeg_finish_decompress(&info) == false) goto err;

        jpeg_destroy_decompress(&info);
        fclose(fHandle);
        delete[] pBuffer;
        return;

    err:
        cout << "(linux) error during: " << szContext << endl;
        exit(-1);
    }

    unique_ptr<Surface> GetSurface()
    {
        return make_unique<Surface>(surf);
    }
};

void app(Image& image); // fwd declaration

int main(int argc, char* argv[])
{
    Image image;
    image.SetFile(argv[1]);
    app(image);
    return 0;
}

#endif // platform implementation

///////////////////
// platform independent stuff

struct Surface : public surf_t
{
    Surface()
    {
        width = height = 0;
        data = NULL;
    }
    Surface(const surf_t& s)
    {
        width = s.width; height = s.height; data = s.data;
    }
    Surface(const uint32_t w, const uint32_t h)
    {
        width = w; height = h;
        data = new uint8_t[width * height](0);
    }
    virtual ~Surface()
    {
        if (data)
            delete[] data;
        data = 0;
    }

    inline uint8_t& At(const int x, const int y) { return data[x + y * width]; }
};

struct Point : public point_t
{
    Point() {}
    Point(const int _x, const int _y)
    {
        x = _x; y = _y;
    }
    Point(const point_t& p)
    {
        x = p.x; y = p.y;
    }
    virtual ~Point() {}

    template< class T >
    T Distance(const point_t& to)
    {
        return sqrt(T(to.x - x) * T(to.x - x) + T(to.y - y) * T(to.y - y));
    }

    // 0 is up, quantized by q
    template< class T >
    T Angle(const point_t& to, const T q)
    {
        auto fa = (atan2(to.y - y, to.x - x) + M_PI) / M_PI / 2;
        return T(fa * q);
    }
};

struct BBox : public bbox_t
{
    BBox() { xmin = xmax = ymin = ymax = 0; }
    BBox(const int xn, const int xm, const int yn, const int ym)
    {
        xmin = xn; xmax = xm; ymin = yn; ymax = ym;
    }
    virtual ~BBox() {}

    void Set(const point_t& p)
    {
        xmin = xmax = p.x;
        ymin = ymax = p.y;
    }

    void Enlarge(const point_t& p)
    {
        xmin = min(xmin, p.x);
        xmax = max(xmax, p.x);
        ymin = min(ymin, p.y);
        ymax = max(ymax, p.y);
    }

    point_t GetCenter()
    {
        return { (xmin + xmax) / 2, (ymin + ymax) / 2 };
    }

};

// https://cppscripts.com/fft-cpp-code
template< class T >
void fft(vector<complex<T>>& x)
{
    size_t N = x.size();
    if (N <= 1) return;

    vector<complex<T>> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < N / 2; ++k) {
        complex<T> t = polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

template< class T >
void NormalizedSpectrum(std::vector<T>& s, const std::vector<T>& v)
{
    vector<complex<T>> x;
    for (auto p : v)
        x.push_back(p);
    fft(x);
    s.clear();
    for (auto p : x)
        s.push_back(norm(p));
}

// with maxing trigger
template< class T >
bool SelfMax(T& m, const T& v)
{
    if (v > m) { m = v; return true; }
    return false;
}

void app(Image& image)
{
    image.Load();

    auto Surf = image.GetSurface();

    cout << image.GetFile() << " width " << Surf.get()->width << " height " << Surf.get()->height << endl;

    // different classifications might be found by tweaking the tunings
    uint8_t tune_black_tol = 128;
    uint8_t tune_adjacency_tol = 4;
    uint8_t tune_pixel_filter_tol = 25;
    uint8_t tune_cloud_filter_factor = 0;
    uint8_t tune_box_size_tol = 25;
    uint8_t tune_circular_frequency_bins = 32; // must be a power of 2
    uint8_t tune_unify_distance = 10;

    // edge/contrast/adjacency filter
    Surface edge(Surf.get()->width, Surf.get()->height);
    for (uint16_t y = 1; y < Surf.get()->height -1; y++)
        for (uint16_t x = 1; x < Surf.get()->width - 1; x++)
        {
            int weight = 0;
            weight += Surf.get()->At(x - 1, y -1) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x + 0, y -1) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x + 1, y -1) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x - 1, y +1) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x + 0, y +1) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x + 1, y +1) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x - 1, y +0) < tune_black_tol ? 1 : 0;
            weight += Surf.get()->At(x + 1, y +0) < tune_black_tol ? 1 : 0;
            edge.At(x, y) = weight > tune_adjacency_tol ? 0xFF : 0x00;
        }

    // point numbering filter
    uint8_t id = 1; // 0 unused
    Surface numb(edge.width, edge.height);
    for (uint16_t y = 1; y < edge.height - 1; y++)
        for (uint16_t x = tune_pixel_filter_tol; x < edge.width - 2; x++)
            if (edge.At(x, y) && !numb.At(x, y))
            {
                uint8_t i = 0, neighbour = 0;
                for (int s = tune_pixel_filter_tol; s > -tune_pixel_filter_tol; s--)
                    if (neighbour = numb.At(x - s, y - 1))
                        if (!i) i = neighbour;
                for (int s = tune_pixel_filter_tol; s > 0; s--)
                    if (neighbour = numb.At(x - s, y))
                        if (!i) i = neighbour;
                if (!i) i = id++;
                numb.At(x, y) = i;
            }

    // clouding / point-grouping
    std::vector< std::vector<Point> > clouds;
    clouds.resize(id);
    for (uint16_t y = 0; y < numb.height; y++)
        for (uint16_t x = 0; x < numb.width; x++)
            if(numb.At(x, y))
                clouds[numb.At(x, y)].emplace_back(point_t({x,y}));

    // small cloud filtering
    // this is kind of arbitrary and perhaps should be eliminated
    if(tune_cloud_filter_factor>0)
    {
        size_t maxcount = 0;
        for (int i = 0; i < id; i++)
            maxcount = max(maxcount,clouds[i].size());
        for (int i = 0; i < id; i++)
            if(clouds[i].size() < max(static_cast<size_t>(10), maxcount / tune_cloud_filter_factor))
                clouds[i].clear();
    }

    // cloud unifications
    // since the point-numbering filter is a TL-DR progressive scanner, shapes have
    // open limbs at the top or left will receive different numberings on those limbs.
    // the process of unification uses a pixel tolerance to join those limbs which have
    // points that are nearby.
    auto fnCheckOtherClouds = [&](int i, int &u)
    {
        for (auto ip : clouds[i])
            for (int j = i + 1; j < id; j++) // j=i+1 for triangular
                for (auto jp : clouds[j])
                    if (ip.Distance<float>(jp) < (float)tune_unify_distance)
                    {
                        clouds[i].insert(std::end(clouds[i]), std::begin(clouds[j]), std::end(clouds[j]));
                        clouds[j].clear();
                        u++;
                        return; // no need to check remaining points in j cloud
                    }

    };
    if(tune_unify_distance>0)
        while(true)
        {
            int unifications = 0;
            for (int i = 0; i < id; i++)
                fnCheckOtherClouds(i, unifications);
            if(!unifications) break;
        }

    // boxing and interior box removal
    std::vector<BBox> bx_limits;
    std::vector<Point> bx_centers;
    bx_limits.resize(id);
    bx_centers.resize(id);
    for(int i = 0; i < id; i++)
        if(clouds[i].size() > 0)
        {
            bx_limits[i].Set(clouds[i].at(0));
            for (point_t p : clouds[i])
                bx_limits[i].Enlarge(p);
            bx_centers[i] = Point(bx_limits[i].GetCenter());
        }
    for (int i = 0; i < id; i++)
        for (int j = 0; j < id; j++)
            if( i != j )
                if (clouds[i].size() > 0 && clouds[j].size() > 0)
                    if (bx_centers[i].Distance<float>(bx_centers[j]) < tune_box_size_tol)
                    {
                        clouds[j].clear();
                        bx_centers[j] = { 0,0 };
                        bx_limits[j] = { 0,0,0,0 };
                    }

    // circularity filter
    std::vector<double> circ_radius;
    circ_radius.resize(id);
    for (int i = 0; i < id; i++)
        if (clouds[i].size() > 0)
        {
            circ_radius[i] = 0;
            for (Point p : clouds[i])
                circ_radius[i] += p.Distance<double>(bx_centers[i]);
            circ_radius[i] /= clouds[i].size();
        }
    std::vector<double> circ_stdev;
    circ_stdev.resize(id);
    for (int i = 0; i < id; i++)
        if (clouds[i].size() > 0)
        {
            circ_stdev[i] = 0;
            for (Point p : clouds[i])
                circ_stdev[i] += pow(circ_radius[i] - p.Distance<double>(bx_centers[i]), 2);
            circ_stdev[i] = pow(circ_stdev[i], 0.5) / clouds[i].size();
        }

    // angular frequency / power filtering
    std::vector< std::vector<double> > angl_power;
    angl_power.resize(id);
    for (int i = 0; i < id; i++)
        if (clouds[i].size() > 0)
        {
            std::vector<double> angl_count;
            angl_count.resize(tune_circular_frequency_bins);
            for (point_t p : clouds[i])
                angl_count.at(bx_centers[i].Angle(p, tune_circular_frequency_bins-1)) += 1;

            std::vector<double> angl_norm_factor;
            angl_norm_factor.resize(tune_circular_frequency_bins);
            for (int j = 0; j < tune_circular_frequency_bins; j++)
                angl_norm_factor.at(j) = angl_count.at(j) / clouds[i].size();

            std::vector<double> angl_excursion;
            angl_excursion.resize(tune_circular_frequency_bins);
            for (Point p : clouds[i])
                SelfMax<double>(angl_excursion.at(bx_centers[i].Angle(p, tune_circular_frequency_bins-1)), p.Distance<double>(bx_centers[i]));

            std::vector<double> angl_excursion_factor;
            angl_excursion_factor.resize(tune_circular_frequency_bins);
            for (int j = 0; j < tune_circular_frequency_bins; j++)
                angl_excursion_factor.at(j) = angl_norm_factor.at(j) * (angl_excursion.at(j) / circ_radius[i]);

            angl_power.resize(tune_circular_frequency_bins);
            NormalizedSpectrum(angl_power[i], angl_excursion_factor);
        }

    // print results
    for (int i = 0; i < id; i++)
        if (clouds[i].size() > 0)
        {
            int sides = 0;
            double side_max = 0;
            for (int j = 2; j < tune_circular_frequency_bins / 2; j++) // nb start at 2 to skip VLF, end at half the F
                if (SelfMax<double>(side_max, angl_power[i].at(j))) sides = j;

            const char* sz[] = { "very likely", "possibly", "likely not" };
            auto p = 1.0 - circ_stdev[i];

            printf(
                "object %d at %4d %4d has %2d sided symmetry and is %s circular\n",
                i, bx_centers[i].x, bx_centers[i].y, sides, sz[ p > .9 ? 0 : (p > .8 ? 1 : 2) ]
            );
        }
}
