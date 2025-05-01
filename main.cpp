// blueseltzer (c) 2025 orthopteroid@gmail.com
// MIT license

// platform headers
#if defined(WIN32)

#include <windows.h>
#include <comdef.h>
#include <combaseapi.h>
#include <wincodec.h>
#include <shlwapi.h>

#elif defined(UNIX)

#endif // platform headers

///////////////////

#include <iostream>
#include <vector>
#include <complex>

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

using namespace std;

struct Surface
{
    uint32_t width, height;
    uint8_t* data;

    Surface()
    {
        width = height = 0;
        data = NULL;
    }
    Surface(const uint32_t w, const uint32_t h)
    {
        width = w; height = h;
        data = NULL;
    }
    virtual ~Surface()
    {
        if (data) free(data);
    }

    void Malloc()
    {
        data = (uint8_t*)malloc(width * height);
        if (data)
            memset(data, 0, width * height);
    }

    inline uint8_t& At(const int x, const int y) { return data[x + y * width]; }
};

struct Point { int x, y; };

template< class T >
T Distance(const Point& a, const Point& b)
{
    return sqrt(T(a.x - b.x) * T(a.x - b.x) + T(a.y - b.y) * T(a.y - b.y));
}

// 0..1, 0 is up
template< class T >
T Angle(const Point& from, const Point& to, const T q)
{
    auto fa = (atan2(to.y - from.y, to.x - from.x) + M_PI) / M_PI / 2;
    return T(fa * q);
}

struct BBox { int xmin, xmax, ymin, ymax; };

void BoxSet(BBox& b, const Point& p)
{
    b.xmin = b.xmax = p.x;
    b.ymin = b.ymax = p.y;
}

void BoxEnlarge(BBox& b, const Point& p)
{
    b.xmin = min(b.xmin, p.x);
    b.xmax = max(b.xmax, p.x);
    b.ymin = min(b.ymin, p.y);
    b.ymax = max(b.ymax, p.y);
}

void BoxCenter(Point& p, const BBox& b)
{
    p.x = (b.xmin + b.xmax) / 2;
    p.y = (b.ymin + b.ymax) / 2;
}

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
    for(auto p: v)
        x.push_back(p);
    fft(x);
    s.clear();
    for (auto p : x)
        s.push_back(norm(p));
}

// with maxing trigger
template< class T >
boolean SelfMax(T& m, const T& v)
{
    if(v > m) { m = v; return true;}
    return false;
}

//////////////////////////

// platform implementation
#if defined(WIN32)

struct Image
{
    const wchar_t* wzFilename;
    const char* szContext;
    HRESULT hr;
    IWICImagingFactory* pIWICFactory;
    IWICBitmapDecoder* pIDecoder;
    IWICBitmapFrameDecode* pIFrame;
    IWICFormatConverter* pIConverter;

    Image()
    {
        wzFilename = NULL;
        szContext = NULL;
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

    void Open()
    {
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

    void GetData(Surface& s)
    {
        szContext = "CreateFormatConverter::GetSize";
        hr = pIConverter->GetSize(&s.width, &s.height);
        if (FAILED(hr)) goto err;

        szContext = "malloc";
        s.Malloc();
        if (!s.data) goto err;

        szContext = "CreateFormatConverter::CopyPixels";
        hr = pIConverter->CopyPixels(NULL, s.width, s.width * s.height, s.data);
        if (FAILED(hr)) goto err;

        return;
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

#elif defined(UNIX)

#endif // platform implementation

////////////////////////////
// platform independent code

void app(Image& image)
{
    image.Open();

    Surface surf;
    image.GetData(surf);

    cout << "image width " << surf.width << " height " << surf.height << endl;

    // different classifications might be found by tweaking the tunings
    uint8_t tune_black_tol = 128;
    uint8_t tune_adjacency_tol = 4;
    uint8_t tune_pixel_filter_tol = 25;
    uint8_t tune_cloud_filter_factor = 10;
    uint8_t tune_box_size_tol = 25;
    uint8_t tune_circular_frequency_bins = 32; // must be a power of 2

    // edge/contrast/adjacency filter
    Surface edge(surf.width, surf.height);
    edge.Malloc();
    for (UINT y = 1; y < surf.height -1; y++)
        for (UINT x = 1; x < surf.width - 1; x++)
        {
            int weight = 0;
            weight += surf.At(x - 1, y -1) < tune_black_tol ? 1 : 0;
            weight += surf.At(x + 0, y -1) < tune_black_tol ? 1 : 0;
            weight += surf.At(x + 1, y -1) < tune_black_tol ? 1 : 0;
            weight += surf.At(x - 1, y +1) < tune_black_tol ? 1 : 0;
            weight += surf.At(x + 0, y +1) < tune_black_tol ? 1 : 0;
            weight += surf.At(x + 1, y +1) < tune_black_tol ? 1 : 0;
            weight += surf.At(x - 1, y +0) < tune_black_tol ? 1 : 0;
            weight += surf.At(x + 1, y +0) < tune_black_tol ? 1 : 0;
            edge.At(x, y) = weight > tune_adjacency_tol ? 0xFF : 0x00;
        }

    // point numbering filter
    uint8_t id = 1; // 0 unused
    Surface numb(edge.width, edge.height);
    numb.Malloc();
    for (UINT y = 1; y < edge.height -1; y++)
        for (UINT x = tune_pixel_filter_tol; x < edge.width - 2; x++)
            if(edge.At(x,y) && !numb.At(x,y))
            {
                uint8_t cat = 0;
                for(int s = tune_pixel_filter_tol; s > -tune_pixel_filter_tol && !cat; s--)
                    if (numb.At(x - s, y - 1)) cat = numb.At(x - s, y - 1);
                for (int s = tune_pixel_filter_tol; s > 0 && !cat; s--)
                    if (numb.At(x - s, y)) cat = numb.At(x - s, y);
                if (!cat)
                    cat = id++;
                numb.At(x, y) = cat;
            }

    // small cloud filtering
    std::vector< std::vector<Point> > clouds;
    clouds.resize(id);
    for (uint16_t y = 0; y < numb.height; y++)
        for (uint16_t x = 0; x < numb.width; x++)
            if(numb.At(x, y))
                clouds[numb.At(x, y)].emplace_back(Point({x,y}));
    size_t maxcount = 0;
    for (int i = 0; i < id; i++)
        maxcount = max(maxcount,clouds[i].size());
    for (int i = 0; i < id; i++)
        if(clouds[i].size() < max(10, maxcount / tune_cloud_filter_factor))
            clouds[i].clear();

    // boxing and interior box removal
    std::vector<BBox> bx_limits;
    std::vector<Point> bx_centers;
    bx_limits.resize(id);
    bx_centers.resize(id);
    for(int i = 0; i < id; i++)
        if(clouds[i].size() > 0)
        {
            BoxSet(bx_limits[i], clouds[i].at(0));
            for (Point p : clouds[i])
                BoxEnlarge(bx_limits[i], p);
            BoxCenter(bx_centers[i], bx_limits[i]);
        }
    for (int i = 0; i < id; i++)
        for (int j = 0; j < id; j++)
            if( i != j )
                if (clouds[i].size() > 0 && clouds[j].size() > 0)
                    if (Distance<float>(bx_centers[i], bx_centers[j]) < tune_box_size_tol)
                    {
                        clouds[j].clear();
                        bx_centers[j] = Point({ 0,0 });
                        bx_limits[j] = BBox({ 0,0,0,0 });
                    }

    // circularity filter
    std::vector<double> circ_radius;
    circ_radius.resize(id);
    for (int i = 0; i < id; i++)
        if (clouds[i].size() > 0)
        {
            circ_radius[i] = 0;
            for (Point p : clouds[i])
                circ_radius[i] += Distance<double>(p, bx_centers[i]);
            circ_radius[i] /= clouds[i].size();
        }
    std::vector<double> circ_stdev;
    circ_stdev.resize(id);
    for (int i = 0; i < id; i++)
        if (clouds[i].size() > 0)
        {
            circ_stdev[i] = 0;
            for (Point p : clouds[i])
                circ_stdev[i] += pow(circ_radius[i] - Distance<double>(p, bx_centers[i]), 2);
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
            for (Point p : clouds[i])
                angl_count.at(Angle(bx_centers[i], p, tune_circular_frequency_bins-1)) += 1;

            std::vector<double> angl_norm_factor;
            angl_norm_factor.resize(tune_circular_frequency_bins);
            for (int j = 0; j < tune_circular_frequency_bins; j++)
                angl_norm_factor.at(j) = angl_count.at(j) / clouds[i].size();

            std::vector<double> angl_excursion;
            angl_excursion.resize(tune_circular_frequency_bins);
            for (Point p : clouds[i])
                SelfMax<double>(angl_excursion.at(Angle(bx_centers[i], p, tune_circular_frequency_bins-1)), Distance<double>(p, bx_centers[i]));

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
