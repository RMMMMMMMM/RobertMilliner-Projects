using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;
using System.ComponentModel.DataAnnotations;
using System;
using System.Linq;
using System.Threading.Tasks;

namespace EarlyBirdnBird.Pages
{
    public class RentalDetailsModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;

        public RentalDetailsModel(ApplicationDbContext dbContext, UserManager<AppUser> userManager)
        {
            _dbContext = dbContext;
            _userManager = userManager;
        }
        public IList<Rental> Rentals { get; set; }
        private int rentalId;

        [BindProperty]
        public InputModel Input { get; set; }

        public class InputModel
        {
            [Required]
            [DataType(DataType.Date)]
            [Display(Name = "Start Date")]
            public DateTime StartTime { get; set; }

            [Required]
            [DataType(DataType.Date)]
            [Display(Name = "End Date")]
            public DateTime EndTime { get; set; }

            [Required]
            [Display(Name = "Guests")]
            public int Guests { get; set; }
        }

		public async Task<IActionResult> OnGetAsync(int? id)
        {
			Rentals = await _dbContext.Rentals
				.Include(r => r.Prices)
                .Include(r => r.Images)
                .Include(r => r.RentalAmenities)
                    .ThenInclude(ra => ra.Amenity)
                .Where(r => r.Id == id)
				.ToListAsync();

			return Page();
		}

        public async Task<IActionResult> OnPostAsync(int? id)
        {

            //NEED ERROR HANDLING!!!!!!!!!!!!!!!!!!!

            var currentUser = await _userManager.GetUserAsync(User);

            rentalId = (int)id;


            // Check if there is an existing reservation within the specified date range and check reservation status
            var existingReservation = await _dbContext.Reservation
                .FirstOrDefaultAsync(r => r.Rental == rentalId && (r.Status == 1 || r.Status == 2) &&
                                           ((Input.StartTime <= r.BookedEnd && Input.EndTime >= r.BookedStart) ||
                                            (Input.EndTime >= r.BookedStart && Input.EndTime <= r.BookedEnd) ||
                                            (Input.StartTime >= r.BookedStart && Input.StartTime <= r.BookedEnd)));

            if (existingReservation != null)
            {
                ModelState.AddModelError(string.Empty, "The rental is unavailable for the specified dates.");
                Rentals = await _dbContext.Rentals
                    .Include(r => r.Prices)
                    .Include(r => r.Images)
                    .Include(r => r.RentalAmenities)
                    .Where(r => r.Id == id)
                    .ToListAsync();
                return Page();
            }



            var Prices = await _dbContext.Prices
                .Where(p => p.RentalId == rentalId).FirstOrDefaultAsync();

            double grandtotal = (double)(Prices.PricePerNight + 10)*(Input.EndTime.Day - Input.StartTime.Day);

            if (currentUser == null)
            {
                return RedirectToPage("/Identity/Account/Login", new { area = "Identity", returnUrl = Url.Page("/RentalDetails", new { id }) });
            }

            Reservation newReservation = new Reservation()
            {
                Rental = rentalId,
                BookedStart = Input.StartTime,
                BookedEnd = Input.EndTime,
                IndividualTotal = Prices.PricePerNight,
                SalesTax = 10, //either we adjust host summary page or this number for taxes, etc.
                GrandTotal = grandtotal,
                Guests = Input.Guests,
                AppUser = currentUser,
                Status = 1 //1 for pending for host approval
            };

            _dbContext.Reservation.Add(newReservation);
            await _dbContext.SaveChangesAsync();

            //return to home
            return RedirectToPage("./Index");
        }
    }
}
