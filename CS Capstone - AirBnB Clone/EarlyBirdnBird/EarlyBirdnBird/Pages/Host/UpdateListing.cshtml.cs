using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using System.ComponentModel.DataAnnotations;

namespace EarlyBirdnBird.Pages.Host
{
    public class UpdateListingModel : PageModel
    {



        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;


        //DB injection
        private readonly UnitOfWork _uow;
        //second injection to make the file path possible for images know where www.root is located
        private readonly IWebHostEnvironment _webHostEnvironment;



        public int rentalID { get; set; }
        public string rentalNAME { get; set; }


        public UpdateListingModel(ApplicationDbContext db, UserManager<AppUser> um, UnitOfWork uow, IWebHostEnvironment webHostEnvironment)
        {
            _dbContext = db;
            _userManager = um;
            _uow = uow;
            _webHostEnvironment = webHostEnvironment;
        }





        [BindProperty]
        public InputModel Input { get; set; }
        //if we bound this as public Rental objRental {get; set; }  i think we woould not need all of these properties/inputs below
        public class InputModel
        {


            // Add image file property
            [Required(ErrorMessage = "Please select an image.")]
            [Display(Name = "Image File")]
            public IFormFile ImageFile { get; set; }

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
            public int? MaxOcupants { get; set; }

            [Display(Name = "Bathrooms")]
            public int? Bathrooms { get; set; }

            [Display(Name = "Beds")]
            public int? Beds { get; set; }



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



        public IActionResult OnGetAsync()
        {


            //i want to prefill the values, due to time, im not doing it



            var temp = Request.Query["rentalID"];

            rentalID = Convert.ToInt32(temp);


            var tempRemtal = (from r in _dbContext.Rentals
                              where r.Id == rentalID
                              select r).FirstOrDefault();

            rentalNAME = tempRemtal.Name;



            ////rental
            //Input.Name = tempRemtal.Name;
            //Input.Address = tempRemtal.Address;
            //Input.City = tempRemtal.City;
            //Input.State = tempRemtal.State;
            //Input.Zip = tempRemtal.Zip;
            //Input.Description = tempRemtal.Description;
            //Input.RentalType = tempRemtal.RentalType;
            //Input.MaxOcupants = tempRemtal.MaxOcupants;
            //Input.Bathrooms = tempRemtal.Bathrooms;
            //Input.Beds = tempRemtal.Beds;
            ////Input.ImageFile = tempRemtal.Images.FirstOrDefault();
            ////amenity
            //Input.HotTub = (from ra in _dbContext.RentalAmenities
            //                where ra.RentalId == rentalID && ra.AmenityId == 1
            //                select ra).Any();
            //Input.BigTV = (from ra in _dbContext.RentalAmenities
            //                where ra.RentalId == rentalID && ra.AmenityId == 2
            //                select ra).Any();
            //Input.NES = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 3
            //               select ra).Any();
            //Input.WiFi = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 4
            //               select ra).Any();
            //Input.SwimmingPool = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 5
            //               select ra).Any();
            //Input.Gym = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 6
            //               select ra).Any();
            //Input.Parking = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 7
            //               select ra).Any();
            //Input.AirConditioning = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 8
            //               select ra).Any();
            //Input.PetFriendly = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 9
            //               select ra).Any();
            //Input.LaundryFacilities = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 10
            //               select ra).Any();
            //Input.Balcony = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 11
            //               select ra).Any();
            //Input.Fireplace = (from ra in _dbContext.RentalAmenities
            //               where ra.RentalId == rentalID && ra.AmenityId == 12
            //               select ra).Any();






            return Page();
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

            var temp = Request.Query["rentalID"];

            rentalID = Convert.ToInt32(temp);


            var tempRemtal = (from r in _dbContext.Rentals
                              where r.Id == rentalID
                              select r).FirstOrDefault();





            tempRemtal.Name = Input.Name;
            tempRemtal.Address = Input.Address;
            tempRemtal.City = Input.City;
            tempRemtal.State = Input.State;
            tempRemtal.Zip = Input.Zip;
            tempRemtal.Description = Input.Description;
            tempRemtal.RentalType = Input.RentalType;
            tempRemtal.MaxOcupants = Input.MaxOcupants;
            tempRemtal.Bathrooms = Input.Bathrooms;
            tempRemtal.Beds = Input.Beds;
            tempRemtal.Published = true;

            List<RentalAmenity> RAlist = new List<RentalAmenity>();
            if (Input.HotTub)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Hot Tub");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.BigTV)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Big TV");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.NES)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "NES");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.WiFi)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "WiFi");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.SwimmingPool)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Swimming Pool");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Gym)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Gym");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Parking)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Parking");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.AirConditioning)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Air Conditioning");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.PetFriendly)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Pet Friendly");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.LaundryFacilities)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Laundry Facilities");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Balcony)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Balcony");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }
            if (Input.Fireplace)
            {
                Amenity tempAmen = _dbContext.Amenities.FirstOrDefault(a => a.Name == "Fireplace");
                RAlist.Add(new RentalAmenity
                {
                    RentalId = tempRemtal.Id,
                    AmenityId = tempAmen.Id,
                    Amenity = tempAmen
                });
            }

            tempRemtal.RentalAmenities = RAlist;
            _dbContext.Update(tempRemtal);
            await _dbContext.SaveChangesAsync();
            ////_dbContext.SaveChanges();

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
                        RentalId = tempRemtal.Id,
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
