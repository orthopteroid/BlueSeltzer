// blueseltzer (c) 2025 orthopteroid@gmail.com
// MIT license

// platform headers
// a little crazy to have here, but it seems MSDEV puts some platforms types into the std headers
#if defined(WIN32)

#include <windows.h>
#include <comdef.h>
#include <combaseapi.h>
#include <wincodec.h>
#include <shlwapi.h>

#elif defined(LINUX)

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
        data = new uint8_t[width * height](0);
    }
    void Realloc(const uint32_t w, const uint32_t h)
    {
        if (data)
            delete[] data;
        width = w; height = h;
        data = new uint8_t[width * height](0);
    }
    virtual ~Surface()
    {
        if (data)
            delete[] data;
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
bool SelfMax(T& m, const T& v)
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

    void Load()
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

    Surface GetData()
    {
        Surface s;
        uint32_t w, h;

        szContext = "CreateFormatConverter::GetSize";
        hr = pIConverter->GetSize(&w, &h);
        if (FAILED(hr)) goto err;

        s.Realloc(w,h);

        szContext = "CreateFormatConverter::CopyPixels";
        hr = pIConverter->CopyPixels(NULL, s.width, s.width * s.height, s.data);
        if (FAILED(hr)) goto err;

        return s;
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

#include <jpeglib.h>
#include <jerror.h>
#include <memory.h>
#include <setjmp.h>

// https://www.tspi.at/2020/03/20/libjpegexample.html#gsc.tab=0
struct Image
{
    struct jpeg_decompress_struct info;
    struct jpeg_error_mgr err;

    const char* szContext;
    uint8_t* pData;

    char* lpFilename;

    Image()
    {
        szContext = 0;
        pData = 0;
    }
    virtual ~Image()
    {
        if (pData)
            delete[] pData;
    }

    void SetFile(char* sz)
    {
        lpFilename = sz;
    }

    void Load()
    {
        FILE* fHandle = NULL;
        uint8_t* pBuffer = 0;

        szContext = "fopen";
        fHandle = fopen(lpFilename, "rb");
        if (fHandle == NULL) goto err;

        info.err = jpeg_std_error(&err);
        jpeg_create_decompress(&info);

        jpeg_stdio_src(&info, fHandle);
        jpeg_read_header(&info, TRUE);

        jpeg_start_decompress(&info);

        pData = new uint8_t[info.output_width * info.output_height];
        pBuffer = new uint8_t[info.output_width * info.num_components];

        while (info.output_scanline < info.output_height)
        {
            uint8_t* ppBuffer[1];
            ppBuffer[0] = pBuffer;
            szContext = "jpeg_read_scanlines";
            JDIMENSION n = jpeg_read_scanlines(&info, ppBuffer, 1);
            if (n!=1) goto err;
            for(int i=0; i<info.output_width; i++)
                pData[info.output_width * (info.output_scanline-1) + i] =
                    (uint8_t)(0.299f * (float)pBuffer[i*3+0] + 0.587f * (float)pBuffer[i*3+1] + 0.114f * (float)pBuffer[i*3+2]);
        }

        jpeg_finish_decompress(&info);
        jpeg_destroy_decompress(&info);
        fclose(fHandle);
        delete[] pBuffer;
        return;

    err:
        cout << "(linux) " << szContext << endl;
        exit(-1);
    }

    Surface GetData()
    {
        Surface s(info.output_width, info.output_height);
        memcpy(s.data, pData, s.width * s.height);
        return s;
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

////////////////////////////
// platform independent code

void app(Image& image)
{
    image.Load();

    Surface surf = image.GetData();

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
    for (uint16_t y = 1; y < surf.height -1; y++)
        for (uint16_t x = 1; x < surf.width - 1; x++)
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
    for (uint16_t y = 1; y < edge.height -1; y++)
        for (uint16_t x = tune_pixel_filter_tol; x < edge.width - 2; x++)
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
        if(clouds[i].size() < max(static_cast<size_t>(10), maxcount / tune_cloud_filter_factor))
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
