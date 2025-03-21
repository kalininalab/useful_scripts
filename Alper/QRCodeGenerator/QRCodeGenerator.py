try:
    import qrcode
    from PIL import Image, ImageDraw, ImageFont
    import sys
except ImportError:
    print("Please install the required libraries using 'pip install -r requirements.txt'")
    exit()

HIPS_COLOUR = "#005aa2"

def qrCodeGenerator(url, savePath, colour=HIPS_COLOUR, back_color="transparent"):
    # URL to encode
    url = url

    # Create QR code
    qr = qrcode.QRCode(
        version=1,
        error_correction=qrcode.constants.ERROR_CORRECT_H,
        box_size=10,
        border=4,
    )
    qr.add_data(url)
    qr.make(fit=True)

    # Generate QR code image
    qr_img = qr.make_image(fill_color=colour, back_color=back_color).convert("RGBA")

    # Save the final image
    if savePath[-4:] != ".png":
        savePath += ".png"
    qr_img.save(f"{savePath}")

if __name__ == "__main__":
    if 3 <= len(sys.argv) <= 5:
        qrCodeGenerator(*sys.argv[1:])
    else:
        print("Usage: python QRCodeGenerator.py <URL> <savePath> [<colour>] [<back_color>]")
        print("Note: <colour> and <back_color> are optional.")
        exit(1)