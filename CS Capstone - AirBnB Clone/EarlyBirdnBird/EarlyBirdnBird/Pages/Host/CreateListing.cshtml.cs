using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Hosting;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.AspNetCore.Mvc.Rendering;
using Microsoft.EntityFrameworkCore;
using SendGrid.Helpers.Mail;
using System;
using System.ComponentModel.DataAnnotations;
using System.Reflection;

namespace EarlyBirdnBird.Pages.Host
{
    public class CreateListingModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;


        //DB injection
        private readonly UnitOfWork _uow;
        //second injection to make the file path possible for images know where www.root is located
        private readonly IWebHostEnvironment _webHostEnvironment;
        
       

        public CreateListingModel(ApplicationDbContext db, UserManager<AppUser> um, UnitOfWork uow, IWebHostEnvironment webHostEnvironment)
        {
            _dbContext = db;
            _userManager = um;
            _uow = uow;
            _webHostEnvironment = webHostEnvironment;
        }
        ///////////// ZACK ADDITION STARTS

        public IActionResult OnGet()
        {
            PriceInput = new PriceInputModel
            {
                PriceDateStart = DateTime.Today,
                PriceDateEnd = DateTime.Today.AddDays(7)
            };

            return Page();
        }

        [BindProperty]
        public PriceInputModel PriceInput { get; set; }
        public class PriceInputModel
        {
            [Required(ErrorMessage = "Please enter the start date.")]
            [Display(Name = "Start Date")]
            public DateTime PriceDateStart { get; set; }

            [Required(ErrorMessage = "Please enter the end date.")]
            [Display(Name = "End Date")]
            public DateTime PriceDateEnd { get; set; }

            [Required(ErrorMessage = "Please enter the price per night.")]
            [Display(Name = "Price Per Night")]
            public double PricePerNight { get; set; }
        }

        ///////////// ZACK ADDITION ENDS
        [BindProperty]
        public InputModel Input { get; set; }
        //if we bound this as public Rental objRental {get; set; }  i think we woould not need all of these properties/inputs below
        public class InputModel
        {


            // Add image file property
            [Required(ErrorMessage = "Please select an image.")]
            [Display(Name = "Image File")]
            public List<IFormFile> ImageFile { get; set; }
            //public IFormFile ImageFile { get; set; }
            

            //rental inputs

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Name")]
            public string Name { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Address")]
            public string Address { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "City")]
            public string City { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "State")]
            public string State { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Zip")]
            public string Zip { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Description")]
            public string Description { get; set; }

            [DataType(DataType.Text)]
            [Display(Name = "RentalType")]
            public string RentalType { get; set; }

            [Display(Name = "MaxOcupants")]
            public int MaxOcupants { get; set; }

            [Display(Name = "Bathrooms")]
            public int Bathrooms { get; set; }

            [Display(Name = "Beds")]
            public int Beds { get; set; }



            //amenity inputs

            [Required]
            [Display(Name = "Hot Tub")]
            public Boolean HotTub { get; set; }

            [Required]
            [Display(Name = "Big TV")]
            public Boolean BigTV { get; set; }

            [Required]
            [Display(Name = "NES")]
            public Boolean NES { get; set; }

            [Required]
            [Display(Name = "WiFi")]
            public Boolean WiFi { get; set; }

            [Required]
            [Display(Name = "Swimming Pool")]
            public Boolean SwimmingPool { get; set; }

            [Required]
            [Display(Name = "Gym")]
            public Boolean Gym { get; set; }

            [Required]
            [Display(Name = "Parking")]
            public Boolean Parking { get; set; }

            [Required]
            [Display(Name = "Air Conditioning")]
            public Boolean AirConditioning { get; set; }

            [Required]
            [Display(Name = "Pet Friendly")]
            public Boolean PetFriendly { get; set; }

            [Required]
            [Display(Name = "Laundry Facilities")]
            public Boolean LaundryFacilities { get; set; }

            [Required]
            [Display(Name = "Balcony")]
            public Boolean Balcony { get; set; }

            [Required]
            [Display(Name = "Fireplace")]
            public Boolean Fireplace { get; set; }


        }





        public async Task<IActionResult> OnPostAsync()
        {

            //NEED ERROR HANDLING!!!!!!!!!!!!!!!!!!!

            var currentUser = await _userManager.GetUserAsync(User);

            if (currentUser == null)
            {
                // User not found, handle the error accordingly
                return NotFound("User not found.");
            }




            Rental newRental = new Rental()
            {
                OwnerName = currentUser,
                Name = Input.Name,
                Address = Input.Address,
                City = Input.City,
                State = Input.State,
                Zip = Input.Zip,
                Description = Input.Description,
                RentalType = Input.RentalType,
                MaxOcupants = Input.MaxOcupants,
                Bathrooms = Input.Bathrooms,
                Beds = Input.Beds,
                Published = true
            };




            _dbContext.Rentals.Add(newRental);
            await _dbContext.SaveChangesAsync();

            ///////////// ZACK ADDITION STARTS
            Price newPrice = new Price
            {
                PriceDateStart = PriceInput.PriceDateStart,
                PriceDateEnd = PriceInput.PriceDateEnd,
                PricePerNight = PriceInput.PricePerNight,
                RentalId = newRental.Id
            };


            _dbContext.Prices.Add(newPrice);
            await _dbContext.SaveChangesAsync();

            //////////// ZACK ADDITION ENDS

            List<RentalAmenity> RAlist = new List<RentalAmenity>();
            if (Input.HotTub)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Hot Tub");
                RAlist.Add(new RentalAmenity {  RentalId = newRental.Id,
                                                AmenityId = tempAmen.Id, 
                                                Amenity = tempAmen });
            }
            if (Input.BigTV)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Big TV");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.NES)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "NES");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.WiFi)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "WiFi");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.SwimmingPool)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Swimming Pool");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Gym)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Gym");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Parking)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Parking");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.AirConditioning)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Air Conditioning");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.PetFriendly)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Pet Friendly");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.LaundryFacilities)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Laundry Facilities");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Balcony)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Balcony");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Fireplace)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Fireplace");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = newRental.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }

            newRental.RentalAmenities = RAlist;
            _dbContext.Update(newRental);
            await _dbContext.SaveChangesAsync();
            //_dbContext.SaveChanges();

            //-----------------------------------------------------------------------------------------------------------
            //determin root path
            string webrootpath = _webHostEnvironment.WebRootPath;
            //retreve files
            var files = HttpContext.Request.Form.Files;
            //image uploaded

            bool isFirstImage = true;
            // Process each uploaded image
            foreach (var imageFile in files)
            {
                if (imageFile.Length > 0)
                {
                    // Create a unique identifier for the image
                    string filename = Guid.NewGuid().ToString();

                    var uploads = Path.Combine(webrootpath, @"images/rentalimg/");

                    // Get and preserve file type
                    var extension = Path.GetExtension(imageFile.FileName);

                    // Create full path of the item to stream
                    var fullPath = Path.Combine(uploads, filename + extension);

                    // Stream binary files to server
                    using var fileStream = System.IO.File.Create(fullPath);
                    await imageFile.CopyToAsync(fileStream);

                    // Associate actual real URL path and save to DB
                    string imagePath = @"/images/rentalimg/" + filename + extension;

                    // Save image information to the database
                    Image newRentalImage = new Image
                    {
                        RentalId = newRental.Id,
                        RentalImage = imagePath,
                        PrimaryImage = isFirstImage // Set PrimaryImage based on isFirstImage flag
                    };

                    _dbContext.Images.Add(newRentalImage);
                    await _dbContext.SaveChangesAsync();

                    // Reset isFirstImage flag after processing the first image
                    isFirstImage = false;
                }
            }

            return RedirectToPage("./Listings");


        }

    }
}
